#ifndef BOX_TYPE_H
#define BOX_TYPE_H
#include<string>
#include<cmath>
#include"uti_func.h"
const double pi = 3.14159265358979;
class box
{
private:
	double a{}, b{}, c{};
	double alpha{}, beta{}, gamma{};
	double v1[3]{}, v2[3]{}, v3[3]{};
	std::string box_type;
public:
	box() :a{}, b{}, c{}, alpha{ 90.0 }, beta{ 90.0 }, gamma{ 90.0 }, box_type{ "unknown" } {};
	template<typename T>
	box(T&& _a, T&& _b, T&& _c, T&& _alpha, T&& _beta, T&& _gamma) :a{ std::forward<T>(_a) }, b{ std::forward<T>(_b) }, c{ std::forward<T>(_c) }
		, alpha{std::forward<T>(_alpha)}, beta{ std::forward<T>(_beta) }, gamma{ std::forward<T>(_gamma) }
	{
		setv();
	}
	box(const box& rhs) :a{ rhs.a }, b{rhs.b }, c{rhs.c}, alpha{rhs.alpha}, beta{rhs.beta}, gamma{rhs.gamma} 
	{
		setv();
	}
	box(box&& rhs) :a{ std::move(rhs.a) }, b{ std::move(rhs.b) }, c{ std::move(rhs.c) }
	, alpha{ std::move(rhs.alpha) }, beta{ std::move(rhs.beta) }, gamma{ std::move(rhs.gamma) } { setv(); }
	box& operator=(const box& box_)
	{
		this->a = box_.a;
		this->b = box_.b;
		this->c = box_.c;
		this->alpha = box_.alpha;
		this->beta = box_.beta;
		this->gamma = box_.gamma;
		this->setv();
		return *this;
	}
	box& operator=(box&& box_){
		this->a = std::move(box_.a);
		this->b = std::move(box_.b);
		this->c = std::move(box_.c);
		this->alpha = std::move(box_.alpha);
		this->beta = std::move(box_.beta);
		this->gamma = std::move(box_.gamma);
		this->setv();
		return *this;
	}
	~box() = default;
	void set_a(double _x){ a = _x; }
	double get_a()const{ return a; }
	void set_b(double _y){ b = _y; }
	double get_b()const{ return b; }
	void set_c(double _z){ c = _z; }
	double get_c()const{ return c; }
	void set_alpha(double _alpha){ alpha = _alpha; }
	double get_alpha()const{ return alpha; }
	void set_beta(double _beta){ beta = _beta; }
	double get_beta()const{ return beta; }
	void set_gamma(double _gamma){ gamma = _gamma; }
	double get_gamma()const{ return gamma; }
	void set_box_type(string _box){ box_type = _box; }
	std::string get_box_type(){ return box_type; }
	void get_v1(double* rhs) { rhs[0] = v1[0];rhs[1] = v1[1];rhs[2] = v1[2]; }
	void get_v2(double* rhs) { rhs[0] = v2[0];rhs[1] = v2[1];rhs[2] = v2[2]; }
	void get_v3(double* rhs) { rhs[0] = v3[0];rhs[1] = v3[1];rhs[2] = v3[2]; }
	void setv() 
	{
		if (alpha == 90.0 && beta == 90.0 && gamma == 90.0 && a == b && b == c)
			box_type = "cubic";
		else if (alpha == 90.0 && beta == 90.0 && gamma == 90.0 && (a != b || a != c))
			box_type = "rectangular";
		else if (alpha == 60.0 && beta == 60.0 && gamma == 90.0)
			box_type = "dodecahedron";
		else if (alpha != 90.0 && beta != 90.0 && gamma != 90.0)
			box_type = "triclinic";
		else
			box_type = "unknown";

		if (box_type == "cubic" || box_type == "rectangular")
		{
			v1[0] = a; v1[1] = 0.0; v1[2] = 0.0;
			v2[0] = 0.0; v2[1] = b; v2[2] = 0.0;
			v3[0] = 0.; v3[1] = 0.0; v3[2] = c;
		}
		else if (box_type == "dodecahedron" || box_type == "triclinic")
		{
			v1[0] = a; v1[1] = 0.0;v1[2] = 0.0;
			v2[0] = b*cos(gamma*pi / 180.0); v2[1] = b*sin(gamma*pi / 180.0);v2[2] = 0.0;
			v3[0] = c*cos(beta*pi / 180.0); v3[1] = c*(cos(alpha*pi / 180.0) - cos(beta*pi / 180.0)*cos(gamma*pi/180.0))/sin(gamma*pi/180.0);
			v3[2] = c*sqrt(1 + 2.0*cos(alpha*pi / 180.0)*cos(beta*pi / 180.0)*cos(gamma*pi / 180.0) - pow(cos(alpha*pi / 180.0), 2.0) - pow(cos(beta*pi / 180.0), 2.0) - pow(cos(gamma*pi / 180.0), 2.0)) / sin(gamma*pi / 180.0);
		}
		else
		{
			v1[0] = 0.0;v1[1] = 0.0;v1[2] = 0.0;
			v2[0] = 0.0;v2[1] = 0.0;v2[2] = 0.0;
			v3[0] = 0;v3[1] = 0.0;v3[2] = 0.0;
		}
	}
	double volume()
	{
		return a * b * c * (1 - square(cos(alpha * pi / 180.)) - square(cos(beta * pi / 180.)) - square(cos(gamma * pi / 180.)))
			+ 2. * (cos(alpha * pi / 180.) * cos(beta * pi / 180.) * cos(gamma * pi / 180.));
	}
	
	int convert_to_v(double _a, double _b, double _c, double* _v1, double* _v2, double* _v3)
	{
		if (box_type == "cubic" || box_type == "rectangular")
		{
			_v1[0] = _a; _v1[1] = 0.0; _v1[2] = 0.0;
			_v2[0] = 0.0; _v2[1] = _b; _v2[2] = 0.0;
			_v3[0] = 0.; _v3[1] = 0.0; _v3[2] = _c;
		}
		else if (box_type == "dodecahedron" || box_type == "triclinic")
		{
			_v1[0] = _a; _v1[1] = 0.0; _v1[2] = 0.0;
			_v2[0] = _b * cos(gamma * pi / 180.0); v2[1] = _b * sin(gamma * pi / 180.0); _v2[2] = 0.0;
			_v3[0] = _c * cos(beta * pi / 180.0); _v3[1] = _c * (cos(alpha * pi / 180.0) - cos(beta * pi / 180.0) * cos(gamma * pi / 180.0)) / sin(gamma * pi / 180.0);
			_v3[2] = _c * sqrt(1 + 2.0 * cos(alpha * pi / 180.0) * cos(beta * pi / 180.0) * cos(gamma * pi / 180.0) - pow(cos(alpha * pi / 180.0), 2.0) - pow(cos(beta * pi / 180.0), 2.0) - pow(cos(gamma * pi / 180.0), 2.0)) / sin(gamma * pi / 180.0);
		}
		else
		{
			_v1[0] = 0.0; _v1[1] = 0.0; _v1[2] = 0.0;
			_v2[0] = 0.0; _v2[1] = 0.0; _v2[2] = 0.0;
			_v3[0] = 0; _v3[1] = 0.0; _v3[2] = 0.0;
		}
		return 0;
	}
};
#endif
