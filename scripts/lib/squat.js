// (c)  oblong industries

'use strict';

/**
 * A simple quaternion library for JavaScript that provides free functions
 * which operate on quadruple arrays.
 *
 * The library assumes quaternions are of the form, given some array
 * of numbers `q`:
 *
 *     q[0] + q[1]*i + q[2]*j + q[3]*k
 *
 * Functions which return a quaternion often have an optional
 * argument, at the end of the argument list, which serves as an "out"
 * parameter.  If the caller passes an object (like an `Array`, or a
 * `Float64Array`) via this argument, the function will set the `'0'`,
 * `'1'`, `'2'`, and `'3'` properties on the object with the computed
 * quaternion's component values.  This can be used to recycle space
 * in a preallocated chunk of memory in an array buffer and avoid
 * allocating space for return values.
 */

window.squat = {};

// internal function that constructs and returns a new, zero
// quaternion
function new_quat() {
  return [0, 0, 0, 0];
}

/**
 * Adds two quaternions.
 */
squat.add = function (q1, q2, out) {
  out = out || new_quat();
  out[0] = q1[0] + q2[0];
  out[1] = q1[1] + q2[1];
  out[2] = q1[2] + q2[2];
  out[3] = q1[3] + q2[3];
  return out;
};

/**
 * Multiplies two quaternions.
 * Note: quaternion multiplication is noncommutative.
 */
squat.mul = function (q1, q2, out) {
  out = out || new_quat();

  var a1 = q1[0], a2 = q2[0],
      b1 = q1[1], b2 = q2[1],
      c1 = q1[2], c2 = q2[2],
      d1 = q1[3], d2 = q2[3];

  out[0] = a1*a2 - b1*b2 - c1*c2 - d1*d2;
  out[1] = a1*b2 + b1*a2 + c1*d2 - d1*c2;
  out[2] = a1*c2 - b1*d2 + c1*a2 + d1*b2;
  out[3] = a1*d2 + b1*c2 - c1*b2 + d1*a2;
  return out;
};

/**
 * Multiplies a quaternion by a scalar.
 */
squat.scale = function (q, x, out) {
  out = out || new_quat();
  var a = q[0], b = q[1], c = q[2], d = q[3];

  out[0] = a*x;
  out[1] = b*x;
  out[2] = c*x;
  out[3] = d*x;
  return out;
};

/**
 * Computes the conjugate of a quaternion.
 */
squat.conjugate = function (q, out) {
  out = out || new_quat();
  var a = q[0], b = q[1], c = q[2], d = q[3];
  out[0] = a;
  out[1] = -b;
  out[2] = -c;
  out[3] = -d;
  return out;
};

/**
 * Computes the inverse, or reciprocal, of a quaternion.
 */
squat.inverse = function (q, out) {
  out = out || new_quat();
  var a = q[0], b = q[1], c = q[2], d = q[3];
  var r = 1 / (a*a + b*b + c*c + d*d);
  out[0] = a*r;
  out[1] = -b*r;
  out[2] = -c*r;
  out[3] = -d*r;
  return out;
};

/**
 * Computes the length of a quaternion: that is, the square root of
 * the product of the quaternion with its conjugate.  Also known as
 * the "norm".
 */
squat.length = function (q) {
  var a = q[0], b = q[1], c = q[2], d = q[3];
  return Math.sqrt(a*a + b*b + c*c + d*d);
};

/**
 * Normalizes a quaternion so its length is equal to 1.  The result of
 * normalizing a zero quaternion is undefined.
 */
squat.normalized = function (q, out) {
  out = out || new_quat();
  var a = q[0], b = q[1], c = q[2], d = q[3];
  var rlen = 1 / squat.length(q);
  return squat.scale(q, rlen, out);
};

/**
 * Provides the real part of the quaternion.
 */
squat.real = function (q) { return q[0]; };

/**
 * Provides the vector part of the quaternion.
 */
squat.vect = function (q) { return [q[1], q[2], q[3]]; };

/**
 * Provides an empty quaternion.
 */
squat.zero = function (out) {
  out = out || new_quat();
  return out;
};

/**
 * Constructs a rotation quaternion from an axis (a normalized
 * "vect3") and an angle (in radians).
 */
squat.from_axis_angle = function (axis, angle, out) {
  out = out || new_quat();
  var x = axis[0], y = axis[1], z = axis[2];
  var r = 1/Math.sqrt(x*x + y*y + z*z);
  var s = Math.sin(angle/2);
  out[0] = Math.cos(angle/2);
  out[1] = s * x * r;
  out[2] = s * y * r;
  out[3] = s * z * r;
  return out;
};

/**
 * Extracts the angle part, in radians, of a rotation quaternion.
 */
squat.angle = function (quat) {
  var a = quat[0];
  if (a < -1.0 || a > 1.0)
    return 0.0;
  var angle = 2 * Math.acos(a);
  if (angle > Math.PI)
    return (angle - 2 * Math.PI);
  return angle;
};

/**
 * Extracts the axis part, as an array of three numbers, of a rotation
 * quaternion.
 */
squat.axis = function (quat) {
  var x = quat[1], y = quat[2], z = quat[3];
  var r = 1/Math.sqrt(x*x + y*y + z*z);
  return [x*r, y*r, z*r];
};

/**
 * Constructs a rotation quaternion from "norm" and "over" vectors.
 */
squat.from_norm_over = function (norm, over, out) {

};

squat.vmul = function(quat, v, out){
    out = out || [0, 0, 0];
 
    var x = v[0],
        y = v[1],
        z = v[2];
 
    var qx = quat[1],
        qy = quat[2],
        qz = quat[3],
        qw = quat[0];
 
    // q*v
    var ix =  qw * x + qy * z - qz * y,
    iy =  qw * y + qz * x - qx * z,
    iz =  qw * z + qx * y - qy * x,
    iw = -qx * x - qy * y - qz * z;
 
    out[0] = ix * qw + iw * -qx + iy * -qz - iz * -qy;
    out[1] = iy * qw + iw * -qy + iz * -qx - ix * -qz;
    out[2] = iz * qw + iw * -qz + ix * -qy - iy * -qx;
 
    return out;
};

squat.crossv = function(a, b, out) {
  out = out || [0, 0, 0];
  const a1 = a[0], a2 = a[1], a3 = a[2];
  const b1 = b[0], b2 = b[1], b3 = b[2];
  out[0] = a2 * b3 - a3 * b2;
  out[1] = a3 * b1 - a1 * b3;
  out[2] = a1 * b2 - a2 * b1;
  return out
};

squat.normv = function(a, out) {
  out = out || [0, 0, 0];
  const length = Math.sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  out[0] = a[0] / length;
  out[1] = a[1] / length;
  out[2] = a[2] / length;
  return out;
}

/**
 * Basis quaternions, for your convenience.
 */
squat.bases = {
  _1: [1, 0, 0, 0],
  i:  [0, 1, 0, 0],
  j:  [0, 0, 1, 0],
  k:  [0, 0, 0, 1]
};