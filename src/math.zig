const std = @import("std");
const math = std.math;

pub const infinity = std.math.inf(f32);
pub const pi = std.math.pi;

pub fn degree_to_radians(deg: f32) f32 {
    return deg * pi / 180.0;
}

pub fn clamp(x: i32, min: i32, max: i32) i32 {
    if (x < min) {
        return min;
    } else if (x > max) {
        return max;
    } else {
        return x;
    }
}

pub const Vector3 = struct {
    x: f32,
    y: f32,
    z: f32,

    pub fn new(x: f32, y: f32, z: f32) Vector3 {
        return Vector3{ .x = x, .y = y, .z = z };
    }

    pub fn random() Vector3 {
        return Vector3.new( random_f32(), random_f32(), random_f32() );
    }

    pub fn random_in_range(min: f32, max: f32) Vector3 {
        return Vector3.new( 
            random_f32_in_range( min, max ), 
            random_f32_in_range( min, max ),
            random_f32_in_range( min, max ) );
    }

    pub fn length(self: *const Vector3) f32 {
        return math.sqrt(self.length_squared());
    }

    pub fn length_squared(self: *const Vector3) f32 {
        return self.x * self.x + self.y * self.y + self.z * self.z;
    }

    pub fn add(self: *Vector3, other: Vector3) void {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }

    pub fn scale(self: *Vector3, s: f32) void {
        self.x *= s;
        self.y *= s;
        self.z *= s;
    }

    pub fn is_near_zero(self: *const Vector3) bool {
        const epsilon = 1e-8;
        return @fabs(self.x) < epsilon and @fabs(self.y) < epsilon and @fabs(self.z) < epsilon;
    }
};

pub fn negate(v: Vector3) Vector3 {
    return Vector3.new(-v.x, -v.y, -v.z);
}

pub fn normalize(v: Vector3) Vector3 {
    const l = v.length();
    return Vector3.new(v.x / l, v.y / l, v.z / l);
}

pub fn reflect(v: Vector3, n: Vector3) Vector3 {
    return sub(v, scale(n, 2 * dot(v, n)));
}

pub fn refract(uv: Vector3, n: Vector3, etai_over_etat: f32) Vector3 {
    const cos_theta = math.min(dot(negate(uv), n), 1);
    const r_out_perp = scale(add(uv, scale(n, cos_theta)), etai_over_etat);
    const r_out_parallel = scale(n, -math.sqrt(@fabs(1 - r_out_perp.length_squared())));
    return add(r_out_perp, r_out_parallel);
}

pub fn lerp(a: Vector3, b: Vector3, t: f32) Vector3 {
    return Vector3.new(t * a.x + (1 - t) * b.x, t * a.y + (1 - t) * b.y, t * a.z + (1 - t) * b.z);
}

pub fn add(a: Vector3, b: Vector3) Vector3 {
    return Vector3.new(a.x + b.x, a.y + b.y, a.z + b.z);
}

pub fn sub(a: Vector3, b: Vector3) Vector3 {
    return Vector3.new(a.x - b.x, a.y - b.y, a.z - b.z);
}

pub fn mul(a: Vector3, b: Vector3) Vector3 {
    return Vector3.new(a.x * b.x, a.y * b.y, a.z * b.z);
}

pub fn scale(a: Vector3, b: f32) Vector3 {
    return Vector3.new(a.x * b, a.y * b, a.z * b);
}

pub fn dot(a: Vector3, b: Vector3) f32 {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

pub fn cross(a: Vector3, b: Vector3) Vector3 {
    return Vector3{
        .x = a.y * b.z - a.z * b.y,
        .y = a.z * b.x - a.x * b.z,
        .z = a.x * b.y - a.y * b.x,
    };
}

pub const Point3 = Vector3;
pub const Color = Vector3;

pub const Ray = struct {
    origin: Point3,
    direction: Vector3,

    pub fn new(origin: Point3, direction: Vector3) Ray {
        return Ray{ .origin = origin, .direction = direction };
    }

    pub fn at(self: *const Ray, t: f32) Vector3 {
        return add(self.origin, scale(self.direction, t));
    }
};

// random numbers:
var g_rng: *std.rand.Random = undefined;

pub fn random_init(rng: *std.rand.Random) void {
    g_rng = rng;
}

pub fn random_f32() f32 {
    return g_rng.float(f32);
}

pub fn random_f32_in_range(min: f32, max: f32) f32 {
    return random_f32() * (max - min) + min;
}

pub fn random_in_cube(min: f32, max: f32) Vector3 {
    return Vector3{ .x = random_f32_in_range(min, max), .y = random_f32_in_range(min, max), .z = random_f32_in_range(min, max) };
}

pub fn random_in_unit_sphere() Vector3 {
    while (true) {
        const p = random_in_cube(-1, 1);

        if (p.length_squared() > 1) {
            continue;
        }

        return p;
    }
}

pub fn random_in_hemisphere(normal: Vector3) Vector3 {
    const in_unit_sphere = random_in_unit_sphere();
    if (dot(in_unit_sphere, normal) > 0.0) {
        return in_unit_sphere;
    } else {
        return negate(in_unit_sphere);
    }
}

pub fn random_in_unit_disk() Vector3 {
    while (true) {
        const p = Vector3.new(random_f32_in_range(-1, 1), random_f32_in_range(-1, 1), 0);
        if (p.length_squared() >= 1) continue;
        return p;
    }
}


