#include <algorithm>
#include <array>
#include <sstream>
#include <utility>
#include <vector>
// for storing coordinates / computing angles
struct vec {
	std::array<double, 3> coords;
	inline vec() : coords{{0, 0, 0}} {}
	inline explicit vec(std::array<double, 3> c) : coords(c) {}
	inline vec(const vec &v) : coords(v.coords) {}
	inline vec(double a, double b, double c) : coords{{a, b, c}} {}
	inline double &xRef() { return coords[0]; }
	inline vec &operator*=(double d) {
		coords[0] *= d;
		coords[1] *= d;
		coords[2] *= d;
		return *this;
	};

	inline vec &operator/=(double d) {
		coords[0] /= d;
		coords[1] /= d;
		coords[2] /= d;
		return *this;
	};

	inline vec &operator+=(const vec &v) {
		coords[0] += v.coords[0];
		coords[1] += v.coords[1];
		coords[2] += v.coords[2];
		return *this;
	}

	inline vec &operator-=(const vec &v) {
		coords[0] -= v.coords[0];
		coords[1] -= v.coords[1];
		coords[2] -= v.coords[2];
		return *this;
	}

	friend inline bool operator==(const vec &v1, const vec &v2);
	friend inline bool operator!=(const vec &v1, const vec &v2);
	friend inline vec operator+(const vec &v1, const vec &v2);
	friend inline vec operator+(const vec &v1, double f);
	friend inline vec operator+(double f, const vec &v1);
	friend inline vec operator-(const vec &v1, const vec &v2);
	friend inline vec operator-(const vec &v1, double f);
	friend inline vec operator-(double f, const vec &v1);
	friend inline vec operator*(double factor, const vec &vector);
	friend inline vec operator*(const vec &vector, double factor);
	friend inline vec operator*(const vec &v1, const vec &v2);
	friend inline vec operator-(const vec &vector);
	friend inline vec operator/(const vec &vector, double divisor);
	friend inline std::ostream &operator<<(std::ostream &out, const vec &v);

	inline double &yRef() { return coords[1]; }
	inline double &zRef() { return coords[2]; }
	inline double x() const { return coords[0]; }
	inline double y() const { return coords[1]; }
	inline double z() const { return coords[2]; }
	inline void setX(const double f) { coords[0] = f; }
	inline void setY(const double f) { coords[1] = f; }
	inline void setZ(const double f) { coords[2] = f; }
	static inline vec zero() { return vec(0.0, 0.0, 0.0); }

	inline vec normalized() const {
		double l = length();
		return l > 0 ? *this / l : zero();
	}
	vec &operator=(const vec &other) {
		if (&other == this) return *this;
		coords = other.coords;
		return *this;
	}
	inline double dot(const vec &v) const {
		return coords[0] * v.coords[0] + coords[1] * v.coords[1] + coords[2] * v.coords[2];
	}
	inline const vec cross(const vec &v) const {
		return vec(coords[1] * v.coords[2] - coords[2] * v.coords[1],
		           coords[2] * v.coords[0] - coords[0] * v.coords[2],
		           coords[0] * v.coords[1] - coords[1] * v.coords[0]);
	}
	double length() const {
		return sqrt(coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2]);
	}
	double sqlength() const {
		return coords[0] * coords[0] + coords[1] * coords[1] + coords[2] * coords[2];
	}
	static double getAngle(const vec &v0, const vec &v1) {
		return acos(std::min<double>(
		    1.0, std::max<double>(-1.0, v0.normalized().dot(v1.normalized()))));
	}
};

inline bool operator==(const vec &v1, const vec &v2) {
	return v1.coords[0] == v2.coords[0] && v1.coords[1] == v2.coords[1] &&
	       v1.coords[2] == v2.coords[2];
}
inline bool operator!=(const vec &v1, const vec &v2) {
	return v1.coords[0] != v2.coords[0] && v1.coords[1] != v2.coords[1] &&
	       v1.coords[2] != v2.coords[2];
}

inline vec operator+(const vec &v1, const vec &v2) {
	return vec(v1.coords[0] + v2.coords[0], v1.coords[1] + v2.coords[1],
	           v1.coords[2] + v2.coords[2]);
}
inline vec operator+(const double f, const vec &v) {
	return vec(v.coords[0] + f, v.coords[1] + f, v.coords[2] + f);
}
inline vec operator+(const vec &v, const double f) {
	return vec(v.coords[0] + f, v.coords[1] + f, v.coords[2] + f);
}

inline vec operator-(const vec &v1, const vec &v2) {
	return vec(v1.coords[0] - v2.coords[0], v1.coords[1] - v2.coords[1],
	           v1.coords[2] - v2.coords[2]);
}
inline vec operator*(const vec &v1, const vec &v2) {
	return vec(v1.coords[0] * v2.coords[0], v1.coords[1] * v2.coords[1],
	           v1.coords[2] * v2.coords[2]);
}

inline vec operator*(const double f, const vec &v) {
	return vec(v.coords[0] * f, v.coords[1] * f, v.coords[2] * f);
}
inline vec operator*(const vec &v, const double f) {
	return vec(v.coords[0] * f, v.coords[1] * f, v.coords[2] * f);
}

inline vec operator-(const double f, const vec &v) {
	return vec(v.coords[0] - f, v.coords[1] - f, v.coords[2] - f);
}
inline vec operator-(const vec &v, const double f) {
	return vec(v.coords[0] - f, v.coords[1] - f, v.coords[2] - f);
}

inline vec operator-(const vec &v) {
	return vec(-v.coords[0], -v.coords[1], -v.coords[2]);
}

inline vec operator/(const vec &v, const double f) {
	return vec(v.coords[0] / f, v.coords[1] / f, v.coords[2] / f);
}

inline std::ostream &operator<<(std::ostream &out, const vec &ve) {
	out << "(" << ve.coords[0] << ", " << ve.coords[1] << ", " << ve.coords[2] << ")";
	return out;
}
