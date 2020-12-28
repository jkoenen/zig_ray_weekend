const std = @import("std");
const math = std.math;
const vec3 = @import("vec3.zig"); // this feels strange.. how are modules done?
const Point = vec3.Point;
const Vec3 = vec3.Vec3;

pub const Ray = struct {
    origin: Point,
    direction: Vec3,

    pub fn new(origin: Point, direction: Vec3) Ray {
        return Ray{ .origin = origin, .direction = direction };
    }

    pub fn at(self: *Ray, t: f32) f32 {
        return add(self.origin, scale(self.direction, t));
    }
};
