const std = @import("std");
usingnamespace @import("math.zig");

const HitInfo = struct {
    position: Point3,
    normal: Vector3,
    t: f32,
    front_face: bool,
    material: *const Material,

    pub fn set_face_normal(self: *HitInfo, ray: Ray, outward_normal: Vector3) void {
        self.front_face = dot(ray.direction, outward_normal) < 0;
        self.normal = if (self.front_face) outward_normal else negate(outward_normal);
    }
};

pub const Material = struct {
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

const DielectricMaterial = struct {
    const Self = @This();
    material: Material,
    index_of_refraction: f32,

    pub fn init(index_of_refraction: f32) Self {
        return Self{
            .material = Material{ .scatterFn = scatter },
            .index_of_refraction = index_of_refraction,
        };
    }

    fn scatter(material: *const Material, ray: Ray, hit_info: HitInfo, attenuation: *Color, scattered_ray: *Ray) bool {
        const self = @fieldParentPtr(Self, "material", material);
        attenuation.* = Color.new(1, 1, 1);

        const refraction_ratio = if (hit_info.front_face) 1.0 / self.index_of_refraction else self.index_of_refraction;

        const unit_direction = normalize(ray.direction);
        const cos_theta = std.math.min(dot(negate(unit_direction), hit_info.normal), 1.0);
        const sin_theta = std.math.sqrt(1.0 - cos_theta * cos_theta);

        const cannot_refract = refraction_ratio * sin_theta > 1.0;
        var direction: Vector3 = undefined;
        if (cannot_refract or reflectance(cos_theta, refraction_ratio) > random_f32()) {
            direction = reflect(unit_direction, hit_info.normal);
        } else {
            direction = refract(unit_direction, hit_info.normal, refraction_ratio);
        }

        scattered_ray.* = Ray.new(hit_info.position, direction);
        return true;
    }

    fn reflectance(cosine: f32, ref_idx: f32) f32 {
        // Use Schlick's approximation for reflectance.
        var r0 = (1.0 - ref_idx) / (1.0 + ref_idx);
        r0 = r0 * r0;
        return r0 + (1.0 - r0) * std.math.pow(f32, (1.0 - cosine), 5.0);
    }
};

const RaySphere = struct {
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

    pub fn intersect(self: *const RaySphere, ray: Ray, t_min: f32, t_max: f32, hit: *HitInfo) bool {
        const oc = sub(ray.origin, self.position);
        const a = ray.direction.length_squared();
        const half_b = dot(oc, ray.direction);
        const c = oc.length_squared() - self.radius * self.radius;

        const discriminant = half_b * half_b - a * c;
        if (discriminant < 0) {
            return false;
        }
        const sqrtd = std.math.sqrt(discriminant);

        var t = (-half_b - sqrtd) / a;
        if (t < t_min or t_max < t) {
            t = (-half_b + sqrtd) / a;

            if (t < t_min or t_max < t) {
                return false;
            }
        }

        const hit_position = ray.at(t);

        const outward_normal = scale(sub(hit_position, self.position), 1.0 / self.radius);
        hit.t = t;
        hit.position = hit_position;
        hit.set_face_normal(ray, outward_normal);
        hit.material = self.material;
        return true;
    }
};

pub const RayScene = struct {
    const SelfType = @This();
    spheres: std.ArrayList(RaySphere),
    lambert_materials: std.ArrayList(LambertMaterial),
    metal_materials: std.ArrayList(MetalMaterial),
    dielectric_materials: std.ArrayList(DielectricMaterial),

    pub fn init(allocator: *std.mem.Allocator, sphere_capacity: usize, material_capacity: usize) !SelfType {
        return RayScene{ 
            .spheres = try std.ArrayList(RaySphere).initCapacity(allocator, sphere_capacity),
            .lambert_materials = try std.ArrayList(LambertMaterial).initCapacity( allocator, material_capacity ),
            .metal_materials = try std.ArrayList(MetalMaterial).initCapacity( allocator, material_capacity ),
            .dielectric_materials = try std.ArrayList(DielectricMaterial).initCapacity( allocator, material_capacity ),
        };
    }

    pub fn deinit(self: *SelfType) void {
        self.spheres.deinit();
        self.lambert_materials.deinit();
        self.metal_materials.deinit();
        self.dielectric_materials.deinit();
    }

    pub fn addLambertMaterial( self: *SelfType, albedo: Color ) *const Material {
        const lambert = self.lambert_materials.addOneAssumeCapacity();
        lambert.* =  LambertMaterial.init( albedo );
        return &lambert.material;
    }

    pub fn addMetalMaterial( self: *SelfType, color: Color, fuzz: f32 ) *const Material {
        const metal = self.metal_materials.addOneAssumeCapacity();
        metal.* = MetalMaterial.init( color, fuzz );
        return &metal.material;        
    }

    pub fn addDielectricMaterial( self: *SelfType, index_of_refraction: f32 ) *const Material {
        const dielectric = self.dielectric_materials.addOneAssumeCapacity();
        dielectric.* = DielectricMaterial.init( index_of_refraction );
        return &dielectric.material;        
    }

    pub fn addSphere(self: *SelfType, position: Point3, radius: f32, material: *const Material) void {
        self.spheres.appendAssumeCapacity(RaySphere.init(position, radius, material));
    }

    fn intersect(self: RayScene, ray: Ray, t_min: f32, t_max: f32, hit_info: *HitInfo) bool {
        var hit_anything = false;
        var t_max_current = t_max;
        for (self.spheres.items) |object| {
            var object_hit: HitInfo = undefined;
            if (object.intersect(ray, t_min, t_max_current, &object_hit)) {
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

    fn skyColor(ray: Ray) Color {
        const unit_direction = normalize(ray.direction);
        const t = 0.5 * (unit_direction.y + 1.0);
        const topColor = Color.new(1, 1, 1);
        return lerp(Color.new(0.5, 0.7, 1), Color.new(1, 1, 1), t);
    }

    pub fn rayColor(self: *const RayScene, ray: Ray, depth: i32) Color {
        if (depth <= 0) {
            return Color{ .x = 0, .y = 0, .z = 0 };
        }

        var hit_info: HitInfo = undefined;
        if (self.intersect(ray, 0.001, infinity, &hit_info)) {
            var scatter_ray: Ray = undefined;
            var attenuation: Color = undefined;

            if (hit_info.material.scatter(ray, hit_info, &attenuation, &scatter_ray)) {
                return mul(attenuation, rayColor(self, scatter_ray, depth - 1));
            } else {
                return Color{ .x = 0, .y = 0, .z = 0 };
            }
        }

        return skyColor(ray);
    }

};

pub const Camera = struct {
    origin: Point3,
    lower_left_corner: Point3,
    horizontal: Vector3,
    vertical: Vector3,
    u: Vector3,
    v: Vector3,
    w: Vector3,
    lens_radius: f32,

    pub fn new(look_from: Point3, look_at: Point3, up: Vector3, vfov: f32, aspect_ratio: f32, aperture: f32, focus_distance: f32) Camera {
        const theta = degree_to_radians(vfov);
        const h = std.math.tan(theta / 2.0);
        const viewport_height: f32 = 2.0 * h;
        const viewport_width = aspect_ratio * viewport_height;

        const w = normalize(sub(look_from, look_at));
        const u = normalize(cross(up, w));
        const v = cross(w, u);

        const origin = look_from;
        const horizontal = scale(u, focus_distance * viewport_width);
        const vertical = scale(v, focus_distance * viewport_height);
        const lower_left_corner = sub(sub(sub(origin, scale(horizontal, 0.5)), scale(vertical, 0.5)), scale(w, focus_distance));

        const lens_radius = aperture / 2.0;

        return Camera{
            .origin = origin,
            .lower_left_corner = lower_left_corner,
            .horizontal = horizontal,
            .vertical = vertical,
            .u = u,
            .v = v,
            .w = w,
            .lens_radius = lens_radius,
        };
    }

    pub fn getRay(self: *const Camera, s: f32, t: f32) Ray {
        const rd = scale(random_in_unit_disk(), self.lens_radius);
        const offset = add(scale(self.u, rd.x), scale(self.v, rd.y));
        const position = add(self.origin, offset);
        return Ray.new(position, Vector3.new(self.lower_left_corner.x + s * self.horizontal.x + t * self.vertical.x - position.x, self.lower_left_corner.y + s * self.horizontal.y + t * self.vertical.y - position.y, self.lower_left_corner.z + s * self.horizontal.z + t * self.vertical.z - position.z));
    }
};

