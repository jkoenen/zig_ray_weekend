const std = @import("std");
const math = std.math;

pub const Vec3 = struct {
    x: f32,
    y: f32,
    z: f32,

    pub fn length(self: *Vec3) f32 {
        return math.sqrt(self.length_squared());
    }
    pub fn length_squared(self: *Vec3) f32 {
        return self.x * self.x + self.y * self.y + self.z * self.z;
    }
};

const Point = Vec3;
const Color = Vec3;
