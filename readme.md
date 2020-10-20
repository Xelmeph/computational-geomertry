F上の2次元平面空間に対するライブラリです．
Fはfloat, double, long doubleであることを期待しています．

何も変更を加えなければ F = double です．



# namespace geometry
## Point
2次元ベクトルを表すデータ構造です．

```using Point = std::complex<double>```

なので基本的には以下のものが使えます．
https://cpprefjp.github.io/reference/complex/complex.html

### I/O
```a b```
の形式での入出力に対応しています．

具体的には，
```cpp
Point x;
std::cin >> x;
```
に ```1 2``` を渡すと x は Point(1,2) と等しくなります．

### 基本演算
* ```F dot(const Point& x, const Point& y)```
* ```F cos(const Point& x, const Point& y)```
* ```F cross(const Point& x, const Point& y)```
* ```F sin(const Point& x, const Point& y)```

x, yのドット積，コサイン，クロス積，サインを返します．

* ```F rad_to_deg(double r)```
* ```F deg_to_rad(double d)```

角度とラジアンの変換を行います．

* ```bool eq(const Point& x, const F& d)```

|x| と d の差がEPS以下であるかを判定します．
## Line / Segment
Line は直線を表し，Segment は線分を表すデータ構造です．

### コンストラクタ
* ``` Line(Point a, Point b)```

二点を通る直線を構成します．

* ``` Segment(Point a, Point b)```

二点を端点とする線分を構成します．

# 執筆途中