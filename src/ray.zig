const std = @import("std");
const math = std.math;
const vec3 = @import("vec3.zig"); // this feels strange.. how are modules done?
const Point = vec3.Point;
const Vec3 = vec3.Vec3;

pub const Ray = struct {
    origin: Point3,
    direction: Vec3,

    pub fn at(t: f32) f32 {
        return add(origin, scale(direction, t));
    }
};
