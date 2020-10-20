#ifndef COMPUTATIONAL_GEOMETRY_HPP
#define COMPUTATIONAL_GEOMETRY_HPP
#include <iostream>
#include <complex>
#include <algorithm>
#include <cmath>

namespace geometry{
    using F = double;
    using Point = std::complex<F>;
    constexpr F PI  = 3.14'159'265'358'979'323'846;
    constexpr F EPS = 1e-9;

    F rad_to_deg(F r){
        return r * 180. / PI;
    }

    F deg_to_rad(F d){
        return d / 180. * PI;
    }

    /// dot product
    F dot(const Point& x, const Point& y){
        return x.real() * y.real() + x.imag() * y.imag();
    }

    /// equal to dot(x/|x|, y/|y|);
    F cos(const Point& x, const Point& y){
        return dot(x/std::abs(x), y/std::abs(y));
    }

    /// cross product
    F cross(const Point& x, const Point& y) {
        return x.real()*y.imag() - x.imag()*y.real();
    }

    /// equal to cross(x/|x|, y/|y|)
    F sin(const Point& x, const Point& y){
        return cross(x/std::abs(x), y/std::abs(y));
    }

    /// input from "x y"
    std::istream& operator>>(std::istream& is, Point& x){
        F a, b;
        is >> a >> b;
        x = Point(a,b);
        return is;
    }

    /// output to "x y"
    std::ostream& operator<<(std::ostream& os, const Point & x){
        return os << x.real() << ' ' << x.imag();
    }

    /// is |x - y| < EPS
    template<class T>
    bool eq(const T&, const F&);

    template<> bool eq(const Point& x, const F& y){
        return std::abs(std::abs(y) - y) < EPS;
    }

    template<> bool eq(const F& x, const F& y){
        return std::abs(y-x) < EPS;
    }

    struct Line{
        Point x, y;
        Point v;

        Line() = default;
        Line(Point a, Point b) :x(a), y(b) {
            v = (y - x);
            v /= std::abs(v);
        }

        Point projection(const Point& z) const {
            return x + dot(z-x, v) * v;
        }

        Point reflection(const Point& z) const {
            return z + (projection(z) - z) * 2.;
        }
    };


    struct Segment : Line{
        Segment() = default;
        Segment(Point a, Point b) : Line(a, b) {}
    };

    /// https://onlinejudge.u-aizu.ac.jp/courses/library/4/CGL/1/CGL_1_C
    /// ccw:counter clockwise, cw:clockwise, olb:online back, olf:online front, os:on segment
    enum Position {ccw = -1, cw = 1, olb = -2, olf = 2, os = 0};

    /// find position of x against p1->p2
    Position where(const Point& p1, const Point& p2, const Point& x){
        auto c1 = cross(p2-p1, x - p1);
        if(c1 > EPS) return Position::ccw;
        if(c1 < -EPS) return Position::cw;
        if(dot(p2-p1, x - p1) < -EPS) return Position::olb;
        if(std::abs(p2-p1) < std::abs(x - p1)) return Position::olf;
        return Position::os;
    }

    Position where(const Segment& L, const Point x){
        return where(L.x, L.y, x);
    }

    /// https://onlinejudge.u-aizu.ac.jp/courses/library/4/CGL/2/CGL_2_B
    bool intersect(const Segment& s1, const Segment& s2){
        int x1 = where(s1.x, s1.y, s2.x);
        int y1 = where(s1.x, s1.y, s2.y);
        int x2 = where(s2.x, s2.y, s1.x);
        int y2 = where(s2.x, s2.y, s1.y);
        return x1 * y1 <= 0 and x2 * y2 <= 0;
    }

    bool parallel(const Line& L, const Line& M){
        return eq(cross(L.v, M.v), 0.);
    }

    bool orthogonal(const Line& L, const Line& M){
        return eq(dot(L.v, M.v), 0.);
    }

    template<class S, class T> F distance(const S&, const T&);
    template<class S, class T> F distance(const S&& s, const T&& t){
        return distance(s, t);
    }

    template<> F distance(const Point& x, const Point& y){
        return std::abs(y-x);
    }

    template<> F distance(const Line& L, const Point& p) {
        return std::abs(p - L.projection(p));
    }

    template<> F distance(const Line& L1, const Line& L2){
        return parallel(L1, L2) ? 0. : distance(L1, L2.x);
    }

    template<> F distance(const Segment& s, const Point& p){
        Point x = s.projection(p);
        if(where(s, x) == Position::os) return distance(p, x);
        return std::min(distance(s.x, p), distance(s.y, p));
    }

    /// https://onlinejudge.u-aizu.ac.jp/courses/library/4/CGL/2/CGL_2_D
    template<> F distance(const Segment& s1, const Segment& s2){
        if(intersect(s1, s2)) return 0;
        return std::min({
            distance(s1, s2.x), distance(s1, s2.y),
            distance(s2, s1.x), distance(s2, s1.y),
        });
    }

    Point cross_point(const Segment& s1, const Segment& s2){
        F d1 = distance(s2, s1.x);
        F d2 = distance(s2, s1.y);
        return s1.x + d1 / (d1+d2) * (s1.y - s1.x);
    }

    class Polygon {
        std::vector<Point> ps;
    public:
        Polygon() = default;
        explicit Polygon(std::vector<Point>& p) :ps(p) {}

        /// sort ps according to counter clockwise
        void sort(){

        }

        Point& operator[](size_t i){
            return ps[i];
        }
    };
}

#endif //COMPUTATIONAL_GEOMETRY_HPP
