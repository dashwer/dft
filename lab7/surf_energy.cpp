#include <array>
#include <cmath>
#include <functional>
#include <stdio.h>

/// PHYSICAL CONSTANTS ///

constexpr double h = 6.62607015 * 1e-34;       ///Постоянная Планка("6.62607015e-34")
constexpr double c = 2.99792458 * 1e8;         ///Скорость света("299792458")
constexpr double e = 1.602176634 * 1e-19;      ///Заряд электрона("1.602176634e-19")
constexpr double m_e = 9.1093837015 * 1e-31;   ///Масса электрона("9.1093837015e-31")
constexpr double pi = 3.14159265358979; 	   ///Число пи("3.14159265358979")

constexpr double a_0 = (h*h * 1e7) / (m_e * 4 * pi*pi * e*e * c*c);        /// Боровский радиус
constexpr double E_h = (m_e * 4*pi*pi * e*e*e*e * c*c*c*c * 1e-14) / (h*h);/// Энергия Хартри 

constexpr double coeff = (E_h*1e3) / (a_0*a_0);

/// SYSTEM PARAMETERS ///

constexpr double n0 = 0.0412;
constexpr int Z = 5;
constexpr double C = 8.82;
constexpr double d = 1.8;
constexpr double r_s = std::pow(4.0*M_PI*n0 / 3.0, -1.0/3.0);
constexpr double r_c = (std::sqrt(2.0) * r_s / 3.0) * std::sqrt(
							0.458 - 2.21 / r_s + 0.9*std::pow(Z, 2.0/3.0)
						  + 0.0071*r_s*r_s / ((1.0 + 0.127*r_s)*(1.0 + 0.127*r_s))
					   );

/// GAUSS METHOD PARAMETERS ///

constexpr int nodes_number = 8;
constexpr double left = 0.0;
constexpr double right = 1.0;


constexpr double weights[nodes_number] = 
{
	0.10122854, 0.22238104,
	0.31370664, 0.36268378,
	0.36268378, 0.31370664,
	0.22238104, 0.10122854
};
	     
constexpr double nodes[nodes_number] = 
{
	-0.96028986, -0.79666648,
	-0.52553242, -0.18343464,
	 0.18343464,  0.52553242,
	 0.79666648,  0.96028986
};

/// HOOKE-JEEVES METHOD PARAMETERS ///

struct Dual
{
	double first, second;
};

constexpr double beta_init = 1.0;
constexpr double delta_init = 0.0;
constexpr Dual h_init = {0.01, 0.01};
constexpr Dual h_lim = {1e-8, 1e-8};

///

using func_t = std::function<double(const double)>;
using func2_t = std::function<double(const Dual)>;
using nodes_container_t = std::array<double, nodes_number>;

double w_kin(const double x);
double w_cul(const double x);
double w_x(const double x);
double w_c(const double x);
double w_xc2(const double x);
double w_kin4(const double x);
double w_xc4(const double x);

struct
{
	double kin, cul, x, c;
	double kin2, xc2;
	double kin4, xc4;
} w;

nodes_container_t nodes_recalculation(const double a, const double b);
double gauss_integral(const func_t& f, const double a, const double b);
double sigma(const Dual& param);
Dual hooke_jeeves(const func2_t& f, const Dual x0);
Dual hooke_jeeves_investigation(const func2_t& f, const Dual x, const Dual h);


int main()
{
	w.kin = gauss_integral(w_kin, left, right);
	w.cul = gauss_integral(w_cul, left, right);
	w.x = gauss_integral(w_x, left, right);
	w.c = gauss_integral(w_c, left, right);
	w.kin2 = n0 * std::log(2.0) / 72.0;
	w.xc2 = gauss_integral(w_xc2, left, right);
	w.kin4 = gauss_integral(w_kin4, left, right);
	w.xc4 = gauss_integral(w_xc4, left, right);

	Dual init_param{beta_init, delta_init};
	Dual param0 = hooke_jeeves(sigma, init_param);
	const double sigma0 = sigma(param0);
	const double beta = param0.first;
	const double delta = param0.second;

	printf("r_c: %5.8f\n", r_c);

	printf("w_kin: %5.8f\n", coeff*w.kin/beta);
	printf("w_cul: %5.8f\n", coeff*w.cul/(beta*beta*beta));
	printf("w_x: %5.8f\n", coeff*w.x/beta);
	printf("w_c: %5.8f\n", coeff*w.c/beta);
	printf("w_kin2: %5.8f\n", coeff*w.kin2*beta);
	printf("w_xc2: %5.8f\n", coeff*w.xc2*beta);
	printf("w_kin4: %5.8f\n", coeff*w.kin4*beta*beta*beta);
	printf("w_xc4: %5.8f\n", coeff*w.xc4*beta*beta*beta);

	printf("Beta: %5.8f\n", beta);
	printf("Delta: %5.8f\n", delta);
	printf("Result: %5.8f\n", coeff*sigma0);
}


double w_kin(const double x)
{
	return 0.3 * std::pow(3*M_PI*M_PI, 2.0/3.0) * std::pow(n0, 5.0/3.0) * (
			   	 	std::pow(1.0 - 0.5*x, 5.0/3.0) + std::pow(0.5*x, 5.0/3.0) - 1.0
		         ) / x;
}


double w_cul(const double x)
{
	return -2.0 * M_PI * n0*n0 * (
			   	 (1.0 - 0.5*x)*(1 - 0.5*x) + 0.25*x*x - (1.0 - 0.5*x)
		    ) / x;
}


double w_x(const double x)
{
	return -0.75 * std::pow(3.0/M_PI, 1.0/3.0) * std::pow(n0, 4.0/3.0) * (
			   	 	 std::pow(1.0 - 0.5*x, 4.0/3.0) + std::pow(0.5*x, 4.0/3.0) - 1.0
		           ) / x;
}


double w_c(const double x)
{
	return -0.056 * std::pow(n0, 4.0/3.0) * (
				
				std::pow(1.0 - 0.5*x, 4.0/3.0) / 
				(0.079 + std::pow(n0, 1.0/3.0) * std::pow(1.0 - 0.5*x, 1.0/3.0))
			  + std::pow(0.5*x, 4.0/3.0) / 
				(0.079 + std::pow(n0, 1.0/3.0) * std::pow(0.5*x, 1.0/3.0))
			  - 1.0 / (0.079 + std::pow(n0, 1.0/3.0))
			   	 
			) / x;
}


double w_xc2(const double x)
{
	return (0.25 * std::pow(n0, 2.0/3.0)) / 
		   (std::pow(M_PI, 5.0/3.0) * std::pow(3.0, 4.0/3.0)) * (
				
				(0.4666 + 0.3735 * std::pow(3.0*M_PI*M_PI*n0 * (1.0 - 0.5*x), -2.0/9.0))
			  * (-0.0085 + 0.3318 * std::pow(3.0*M_PI*M_PI*n0 * (1.0 - 0.5*x), 1.0/15.0)) 
			  * (-0.0085 + 0.3318 * std::pow(3.0*M_PI*M_PI*n0 * (1.0 - 0.5*x), 1.0/15.0))
			  / std::pow(1.0 - 0.5*x, 4.0/3.0)

			  + (0.4666 + 0.3735 * std::pow(3.0*M_PI*M_PI*n0 * 0.5*x, -2.0/9.0))
			  * (-0.0085 + 0.3318 * std::pow(3.0*M_PI*M_PI*n0 * 0.5*x, 1.0/15.0)) 
			  * (-0.0085 + 0.3318 * std::pow(3.0*M_PI*M_PI*n0 * 0.5*x, 1.0/15.0))
			  / std::pow(0.5*x, 4.0/3.0)
			   	 
			) * x;
}


double w_kin4(const double x)
{
	using std::pow;

	return (1.336 * pow(n0, 1.0/3.0) / (2060 * pow(3.0 * M_PI*M_PI, 2.0/3.0)))
	     * (
	     	pow(1 - 0.5*x, -11.0/3.0) * ((1 - 0.5*x)*(1 + x/16.0) + x*x/3.0)
	      + 29.0 * x*x * pow(0.5*x, -11.0/3.0) / 96.0
	     ) * x; 
}


double w_xc4(const double x)
{
	using std::pow, std::exp;

	return 7.35 * pow(10.0, -5.0) * (
			x * exp(-0.2986 * pow(n0 * (1 - 0.5*x), -0.26)) / ((1 - 0.5*x)*(1 - 0.5*x))
		  + 4.0 * exp(-0.2986 * pow(n0 * 0.5*x, -0.26)) / x
		);
}


nodes_container_t nodes_recalculation(
	const double a, 
	const double b
)
{
	nodes_container_t output_nodes;

	for (int i = 0; i < nodes_number; ++i) {
		output_nodes[i] = a + 0.5 * (b - a) * (1.0 + nodes[i]);
	}

	return output_nodes;
}


double gauss_integral(
	const func_t& f, 
	const double a, 
	const double b
)
{
	nodes_container_t recalculated_nodes = nodes_recalculation(a, b);
	double I = 0.0;

	for (int i = 0; i < nodes_number; ++i) {
		I += weights[i] * f(recalculated_nodes[i]);
	}
	I *= (b - a) / 2.0;

	return I;
}


double sigma(const Dual& param)
{
	using std::sqrt;
	using std::exp;
	using std::cosh;

	const double beta = param.first;
	const double delta = param.second;

	return w.kin / beta 
		 + w.cul / (beta*beta*beta)
		 + w.x / beta + w.c / beta
		 + w.kin2 * beta + w.xc2 * beta
		 + sqrt(3.0)*Z*Z / (C*C*C) * exp(-4.0*M_PI*d / (sqrt(3.0)*C)) * exp(8.0*M_PI*delta/(sqrt(3.0)*C))
		 + 2*M_PI*n0*n0 * (
		       1.0 - cosh(beta*r_c)*beta*d*exp(-0.5*beta*d) / (1.0-exp(-beta*d))
		   ) / (beta*beta*beta)
		 + 2.0 * M_PI * n0*n0 * d * (
		 	delta*delta + (1 - exp(beta*delta)) * exp(-0.5*beta*d) * cosh(beta*r_c) / (beta*beta)
		 ) 
		 + w.kin4 * (beta*beta*beta)
		 + w.xc4 * (beta*beta*beta);
}


bool operator>(const Dual& lhs, const Dual& rhs)
{
	return (lhs.first > rhs.first) && (lhs.second > rhs.second);
}

Dual operator/=(Dual& lhs, const double rhs)
{
	lhs.first /= rhs;
	lhs.second /= rhs;

	return lhs;
}

Dual operator+(const Dual& lhs, const Dual& rhs)
{
	return Dual{lhs.first + rhs.first, lhs.second + rhs.second};
}

Dual operator-(const Dual& lhs, const Dual& rhs)
{
	return Dual{lhs.first - rhs.first, lhs.second - rhs.second};
}

Dual operator*(const double lhs, const Dual& rhs)
{
	return Dual{lhs * rhs.first, lhs * rhs.second};
}

Dual hooke_jeeves(
	const func2_t& f, 
	const Dual x0
)
{
	Dual h = h_init;
	Dual x_prev = x0, x_next;

	while (h > h_lim) {
		x_next = hooke_jeeves_investigation(f, x_prev, h);

		if (f(x_next) >= f(x_prev)) {	
			h /= 10;
			continue;
		}

		Dual new_x_prev = x_prev + 2 * (x_next - x_prev);
		Dual new_x_next = hooke_jeeves_investigation(f, new_x_prev, h);
		while (f(new_x_next) < f(x_next)) {
			x_prev = new_x_prev;
			x_next = new_x_next;

			new_x_prev = x_prev + 2 * (x_next - x_prev);
			new_x_next = hooke_jeeves_investigation(f, new_x_prev, h);
		}

		h /= 10;
	}

	return x_next;
}


Dual hooke_jeeves_investigation(
	const func2_t& f, 
	const Dual x,
	const Dual h
)
{
	Dual new_x;

	if (f(x + h) < f(x)) {
		new_x = x + h;
	} else if (f(x - h) < f(x)) {
		new_x = x - h;
	} else {
		new_x = x;
	}

	return new_x;
}