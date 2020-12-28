const std = @import("std");
const math = std.math;

pub const Vec3 = struct {
    x: f32,
    y: f32,
    z: f32,

    pub fn new(x: f32, y: f32, z: f32) Vec3 {
        return Vec3{ .x = x, .y = y, .z = z };
    }

    pub fn length(self: *const Vec3) f32 {
        return math.sqrt(self.length_squared());
    }
    pub fn length_squared(self: *const Vec3) f32 {
        return self.x * self.x + self.y * self.y + self.z * self.z;
    }
};

pub fn normalize(v: Vec3) Vec3 {
    const l = v.length();
    return Vec3.new(v.x / l, v.y / l, v.z / l);
}

pub fn lerp(a: Vec3, b: Vec3, t: f32) Vec3 {
    return Vec3.new(t * a.x + (1 - t) * b.x, t * a.y + (1 - t) * b.y, t * a.z + (1 - t) * b.z);
}

pub fn add(a: Vec3, b: Vec3) Vec3 {
    return Vec3.new(a.x + b.x, a.y + b.y, a.z + b.z);
}

pub fn sub(a: Vec3, b: Vec3) Vec3 {
    return Vec3.new(a.x - b.x, a.y - b.y, a.z - b.z);
}

pub fn scale(a: Vec3, b: f32) Vec3 {
    return Vec3.new(a.x * b, a.y * b, a.z * b);
}

pub fn dot(a: Vec3, b: Vec3) f32 {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

pub const Point = Vec3;
pub const Color = Vec3;
