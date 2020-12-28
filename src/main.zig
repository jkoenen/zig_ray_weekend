// Ray Tracing in One Weekend:
// https://raytracing.github.io/books/RayTracingInOneWeekend.html#overview
// lets test zig..

// open questions:
// - how to switch between debug/release build configs?
//      -> debug is default, -O ReleaseFast switches to optimized build

const std = @import("std");
const math = std.math;

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

    pub fn set_face_normal(self: *HitInfo, r: Ray, outward_normal: Vector3) void {
        self.front_face = dot(r.direction, outward_normal) < 0;
        self.normal = if (self.front_face) outward_normal else negate(outward_normal);
    }
};

const HitSphere = struct {
    const SelfType = @This();
    position: Point3,
    radius: f32,

    pub fn new(position: Point3, radius: f32) HitSphere {
        var result: HitSphere = undefined;
        result.init(position, radius);
        return result;
    }

    pub fn init(self: *SelfType, position: Point3, radius: f32) void {
        self.position = position;
        self.radius = radius;
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

    pub fn addSphere(self: *SelfType, position: Point3, radius: f32) !void {
        const sphere = try self.spheres.addOne();
        sphere.init(position, radius);
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

fn rayColor(hit_object: *const HitObjectList, r: Ray) Color {
    var hit_info: HitInfo = undefined;
    if (hit_object.intersect(r, 0, infinity, &hit_info)) {
        return Color.new(0.5 * (hit_info.normal.x + 1), 0.5 * (hit_info.normal.y + 1), 0.5 * (hit_info.normal.z + 1));
    }

    return skyColor(r);
}

fn writePpmColor(out: anytype, color: Color) !void {
    const ir = @floatToInt(i32, 255.0 * color.x + 0.5);
    const ig = @floatToInt(i32, 255.0 * color.y + 0.5);
    const ib = @floatToInt(i32, 255.0 * color.z + 0.5);

    try std.fmt.format(out, "{} {} {}\n", .{ ir, ig, ib });
}

pub fn main() anyerror!void {
    const stdout = std.io.getStdOut().writer();
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    var allocator = &gpa.allocator;

    const blub = try allocator.create(HitSphere);

    var world = HitObjectList.new(allocator);
    defer world.deinit();

    try world.addSphere(Point3.new(0, 0, -1), 0.5);
    try world.addSphere(Point3.new(0, -100.5, -1), 100);

    const hit_object = &world;

    const aspect_ratio: f32 = 16.0 / 9.0;
    const image_width: i32 = 256;
    const image_height: i32 = 256 / aspect_ratio;

    const viewport_height: f32 = 2.0;
    const viewport_width: f32 = aspect_ratio * viewport_height;
    const focal_length: f32 = 1;

    const origin = Vector3.new(0, 0, 0);
    const horizontal = Vector3.new(viewport_width, 0, 0);
    const vertical = Vector3.new(0, viewport_height, 0);
    const lower_left_corner = Vector3.new(origin.x - horizontal.x / 2 - vertical.x / 2, origin.y - horizontal.y / 2 - vertical.y / 2, origin.z - focal_length);

    const ppm_file = try std.fs.cwd().createFile("test.ppm", .{});
    defer ppm_file.close();

    const ppm_writer = ppm_file.writer();

    try ppm_writer.print("P3\n{} {}\n255\n", .{ image_width, image_height });

    // hmm.. how do for-loops work in zig? -> using while loops for now ;)
    var y: i32 = image_height - 1;
    while (y >= 0) {
        try stdout.print("\rScanlines remaining: {}", .{y});
        var x: i32 = 0;
        while (x < image_width) {
            const u = @intToFloat(f32, x) / @intToFloat(f32, image_width - 1);
            const v = @intToFloat(f32, y) / @intToFloat(f32, image_height - 1);

            const r = Ray.new(origin, Vector3.new(lower_left_corner.x + u * horizontal.x + v * vertical.x - origin.x, lower_left_corner.y + u * horizontal.y + v * vertical.y - origin.y, lower_left_corner.z + u * horizontal.z + v * vertical.z - origin.z));

            const pixel_color = rayColor(hit_object, r);
            try writePpmColor(ppm_writer, pixel_color);

            x += 1;
        }
        y -= 1;
    }

    try stdout.print("\nDone\n", .{});
}
