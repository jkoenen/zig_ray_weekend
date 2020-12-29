// Ray Tracing in One Weekend:
// https://raytracing.github.io/books/RayTracingInOneWeekend.html#overview
// lets test zig..

// open questions:
// - how to switch between debug/release build configs?
//      -> debug is default, -O ReleaseFast switches to optimized build

const std = @import("std");
const math = std.math;
const rand = std.rand;
const Random = rand.Random;

const ray_math = @import("math.zig");
usingnamespace ray_math;

fn createNormalizedPosition(x: i32, y: i32, width: i32, height: i32, z: f32) Vector3 {
    // is there a nicer / simpler way for this?
    const xf = @intToFloat(f32, x) / @intToFloat(f32, width - 1);
    const yf = @intToFloat(f32, y) / @intToFloat(f32, height - 1);
    const zf = z;

    return .{ .x = xf, .y = yf, .z = zf };
}

fn skyColor(r: Ray) Color {
    const unit_direction = normalize(r.direction);
    const t = 0.5 * (unit_direction.y + 1.0);
    const topColor = Color.new(1, 1, 1);
    return lerp(Color.new(0.5, 0.7, 1), Color.new(1, 1, 1), t);
}

const HitInfo = struct {
    position: Point3,
    normal: Vector3,
    t: f32,
    front_face: bool,
    material: *const Material,

    pub fn set_face_normal(self: *HitInfo, r: Ray, outward_normal: Vector3) void {
        self.front_face = dot(r.direction, outward_normal) < 0;
        self.normal = if (self.front_face) outward_normal else negate(outward_normal);
    }
};

const Material = struct {
    scatterFn: fn (material: *const Material, ray: Ray, hit_info: HitInfo, attenuation: *Color, scattered_ray: *Ray) bool,

    pub fn scatter(self: *const Material, ray: Ray, hit_info: HitInfo, attenuation: *Color, scattered_ray: *Ray) bool {
        return self.scatterFn(self, ray, hit_info, attenuation, scattered_ray);
    }
};

const LambertMaterial = struct {
    const Self = @This();
    material: Material,
    albedo: Color,

    pub fn init(albedo: Color) Self {
        return Self{
            .material = Material{ .scatterFn = scatter },
            .albedo = albedo,
        };
    }

    fn scatter(material: *const Material, ray: Ray, hit_info: HitInfo, attenuation: *Color, scattered_ray: *Ray) bool {
        const self = @fieldParentPtr(Self, "material", material);
        var scatter_direction = random_in_hemisphere(hit_info.normal);

        // catch degenerate scatter direction:
        if (scatter_direction.is_near_zero()) {
            scatter_direction = hit_info.normal;
        }

        scattered_ray.* = Ray.new(hit_info.position, scatter_direction);
        attenuation.* = self.albedo;
        return true;
    }
};

const MetalMaterial = struct {
    const Self = @This();
    material: Material,
    color: Color,
    fuzz: f32,

    pub fn init(color: Color, fuzz: f32) Self {
        return Self{
            .material = Material{ .scatterFn = scatter },
            .color = color,
            .fuzz = fuzz,
        };
    }

    fn scatter(material: *const Material, ray: Ray, hit_info: HitInfo, attenuation: *Color, scattered_ray: *Ray) bool {
        const self = @fieldParentPtr(Self, "material", material);
        var reflected = reflect(normalize(ray.direction), hit_info.normal);
        scattered_ray.* = Ray.new(hit_info.position, add(reflected, scale(random_in_unit_sphere(), self.fuzz)));
        attenuation.* = self.color;
        return true;
    }
};

const HitSphere = struct {
    const SelfType = @This();
    position: Point3,
    radius: f32,
    material: *const Material,

    pub fn init(position: Point3, radius: f32, material: *const Material) SelfType {
        return SelfType{
            .position = position,
            .radius = radius,
            .material = material,
        };
    }

    pub fn intersect(self: *const HitSphere, r: Ray, t_min: f32, t_max: f32, hit: *HitInfo) bool {
        const oc = sub(r.origin, self.position);
        const a = r.direction.length_squared();
        const half_b = dot(oc, r.direction);
        const c = oc.length_squared() - self.radius * self.radius;
        const discriminant = half_b * half_b - a * c;
        if (discriminant < 0) {
            return false;
        }

        const sqrtd = math.sqrt(discriminant);
        var t = (-half_b - sqrtd) / a;

        if (t < t_min or t_max < t) {
            t = (-half_b + sqrtd) / a;

            if (t < t_min or t_max < t) {
                return false;
            }
        }

        const hit_position = r.at(t);

        const outward_normal = normalize(sub(hit_position, self.position));
        hit.t = t;
        hit.position = hit_position;
        hit.set_face_normal(r, outward_normal);
        hit.material = self.material;
        return true;
    }
};

const HitObjectList = struct {
    const SelfType = @This();
    spheres: std.ArrayList(HitSphere),

    pub fn new(allocator: *std.mem.Allocator) SelfType {
        return HitObjectList{ .spheres = std.ArrayList(HitSphere).init(allocator) };
    }

    pub fn deinit(self: *SelfType) void {
        self.spheres.deinit();
    }

    pub fn addSphere(self: *SelfType, position: Point3, radius: f32, material: *const Material) !void {
        try self.spheres.append(HitSphere.init(position, radius, material));
    }

    pub fn intersect(self: HitObjectList, r: Ray, t_min: f32, t_max: f32, hit_info: *HitInfo) bool {
        var hit_anything = false;
        var t_max_current = t_max;
        for (self.spheres.items) |object| {
            var object_hit: HitInfo = undefined;
            if (object.intersect(r, t_min, t_max_current, &object_hit)) {
                if (hit_anything) {
                    std.debug.assert(object_hit.t <= t_max_current);
                }
                hit_anything = true;
                t_max_current = object_hit.t;
                hit_info.* = object_hit;
            }
        }
        return hit_anything;
    }
};

const Camera = struct {
    origin: Point3,
    lower_left_corner: Point3,
    horizontal: Vector3,
    vertical: Vector3,

    pub fn new(width: i32, height: i32) Camera {
        const aspect_ratio: f32 = @intToFloat(f32, width) / @intToFloat(f32, height);

        const viewport_height: f32 = 2.0;
        const viewport_width: f32 = aspect_ratio * viewport_height;
        const focal_length: f32 = 1;

        //const origin = Vector3.new(0, 0.5, 0.5);
        const origin = Vector3.new(0, 0.0, 0.0);
        const horizontal = Vector3.new(viewport_width, 0, 0);
        const vertical = Vector3.new(0, viewport_height, 0);
        const lower_left_corner = Vector3.new(origin.x - horizontal.x / 2 - vertical.x / 2, origin.y - horizontal.y / 2 - vertical.y / 2, origin.z - focal_length);

        return Camera{
            .origin = origin,
            .lower_left_corner = lower_left_corner,
            .horizontal = horizontal,
            .vertical = vertical,
        };
    }

    pub fn getRay(self: *const Camera, u: f32, v: f32) Ray {
        return Ray.new(self.origin, Vector3.new(self.lower_left_corner.x + u * self.horizontal.x + v * self.vertical.x - self.origin.x, self.lower_left_corner.y + u * self.horizontal.y + v * self.vertical.y - self.origin.y, self.lower_left_corner.z + u * self.horizontal.z + v * self.vertical.z - self.origin.z));
    }
};

fn clamp(x: i32, min: i32, max: i32) i32 {
    if (x < min) {
        return min;
    } else if (x > max) {
        return max;
    } else {
        return x;
    }
}

fn rayColor(world: *const HitObjectList, r: Ray, depth: i32) Color {
    if (depth <= 0) {
        return Color{ .x = 0, .y = 0, .z = 0 };
    }

    var hit_info: HitInfo = undefined;
    if (world.intersect(r, 0.001, infinity, &hit_info)) {
        var scatter_ray: Ray = undefined;
        var attenuation: Color = undefined;

        if (hit_info.material.scatter(r, hit_info, &attenuation, &scatter_ray)) {
            //const target = add(add(hit_info.position, hit_info.normal), random_in_unit_sphere(g_rng));
            return mul(attenuation, rayColor(world, scatter_ray, depth - 1));
        } else {
            return Color{ .x = 0, .y = 0, .z = 0 };
        }
    }

    return skyColor(r);
}

fn writePpmColor(out: anytype, color: Color) !void {
    const g_r = math.sqrt(color.x);
    const g_g = math.sqrt(color.y);
    const g_b = math.sqrt(color.z);
    const ir = clamp(@floatToInt(i32, 255.0 * g_r + 0.5), 0, 255);
    const ig = clamp(@floatToInt(i32, 255.0 * g_g + 0.5), 0, 255);
    const ib = clamp(@floatToInt(i32, 255.0 * g_b + 0.5), 0, 255);

    try std.fmt.format(out, "{} {} {}\n", .{ ir, ig, ib });
}

pub fn main() anyerror!void {
    const stdout = std.io.getStdOut().writer();

    // allocator:
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    var allocator = &gpa.allocator;

    // init random numbers:
    const time_now = std.time.timestamp();
    const seed_i64 = if (time_now >= 0) time_now else -time_now;
    const seed_u64 = @intCast(u64, seed_i64);
    var rng = rand.DefaultPrng.init(seed_u64);
    random_init(&rng.random);

    // create the scene/world:
    var world = HitObjectList.new(allocator);
    defer world.deinit();

    if (false) {
        const matGreen = LambertMaterial.init(Color.new(0.2, 0.6, 0.1));
        const matRed = LambertMaterial.init(Color.new(0.7, 0.2, 0.1));
        const matMetal = MetalMaterial.init(Color.new(0.7, 0.7, 0.7));

        try world.addSphere(Point3.new(0, 0, -1), 0.5, &matRed.material);
        try world.addSphere(Point3.new(0, 0.7, -1), 0.3, &matMetal.material);
        try world.addSphere(Point3.new(0, 1.1, -1), 0.2, &matRed.material);
        try world.addSphere(Point3.new(0, -100.5, -1), 100, &matGreen.material);
    } else {
        const material_ground = LambertMaterial.init(Color.new(0.8, 0.8, 0.0));
        const material_center = LambertMaterial.init(Color.new(0.7, 0.3, 0.3));
        const material_left = MetalMaterial.init(Color.new(0.8, 0.8, 0.8), 0.3);
        const material_right = MetalMaterial.init(Color.new(0.8, 0.6, 0.2), 1.0);

        try world.addSphere(Point3.new(0.0, -100.5, -1.0), 100.0, &material_ground.material);
        try world.addSphere(Point3.new(0.0, 0.0, -1.0), 0.5, &material_center.material);
        try world.addSphere(Point3.new(-1.0, 0.0, -1.0), 0.5, &material_left.material);
        try world.addSphere(Point3.new(1.0, 0.0, -1.0), 0.5, &material_right.material);
    }

    const resolution_scale: i32 = 4;
    const samples_per_pixel: i32 = 100;
    const max_depth: i32 = 50;

    const image_width = 1280 / resolution_scale;
    const image_height = 720 / resolution_scale;

    // camera settings:
    const camera = Camera.new(image_width, image_height);

    // open the output file:
    const ppm_file = try std.fs.cwd().createFile("test.ppm", .{});
    defer ppm_file.close();

    const ppm_writer = ppm_file.writer();

    try ppm_writer.print("P3\n{} {}\n255\n", .{ image_width, image_height });

    // hmm.. how do for-loops work in zig? -> using while loops for now ;)
    var y: i32 = image_height - 1;
    while (y >= 0) {
        try stdout.print("\rScanlines remaining: {} ", .{y});
        var x: i32 = 0;
        while (x < image_width) {
            var pixel_color = Color{ .x = 0, .y = 0, .z = 0 };

            var sample_index: i32 = 0;
            while (sample_index < samples_per_pixel) {
                const u = (@intToFloat(f32, x) + rng.random.float(f32)) / @intToFloat(f32, image_width - 1);
                const v = (@intToFloat(f32, y) + rng.random.float(f32)) / @intToFloat(f32, image_height - 1);

                const ray = camera.getRay(u, v);

                const sample_color = rayColor(&world, ray, max_depth);

                pixel_color.add(sample_color);
                sample_index += 1;
            }
            pixel_color.scale(1.0 / @intToFloat(f32, samples_per_pixel));
            try writePpmColor(ppm_writer, pixel_color);

            x += 1;
        }
        y -= 1;
    }

    try stdout.print("\nDone\n", .{});
}
