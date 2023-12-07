#include <array>
#include <cmath>
#include <functional>
#include <stdio.h>

using w_t = std::function<double(const double)>;
using sigma_t = std::function<double(const double, const double)>;

/// PHYSICAL CONSTANTS ///

constexpr double h = 6.62607015 * 1e-34;      
constexpr double c = 2.99792458 * 1e8;         
constexpr double e = 1.602176634 * 1e-19;      
constexpr double m_e = 9.1093837015 * 1e-31;   
constexpr double pi = 3.14159265358979; 	   

constexpr double a_0 = (h*h * 1e7) / (m_e * 4 * pi*pi * e*e * c*c);         
constexpr double E_h = (m_e * 4*pi*pi * e*e*e*e * c*c*c*c * 1e-14) / (h*h);

constexpr double to_mJ = (E_h*1e3) / (a_0*a_0);
constexpr double to_eV = E_h / e;

/// SYSTEM PARAMETERS ///

constexpr double n0 = 0.0194;
constexpr int Z = 4;
constexpr double C = 6.589;
constexpr double d = 4.66;
constexpr double r_s = std::pow(4.0*M_PI*n0 / 3.0, -1.0/3.0);
constexpr double r_c = (std::sqrt(2.0) * r_s / 3.0) * std::sqrt(
							0.458 - 2.21 / r_s + 0.9*std::pow(Z, 2.0/3.0)
						  + 0.0071*r_s*r_s / ((1.0 + 0.127*r_s)*(1.0 + 0.127*r_s))
					   );
/// GAUSS METHOD ///

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

using nodes_container_t = std::array<double, nodes_number>;

nodes_container_t nodes_recalculation(const double a, const double b);
double gauss_integral(const w_t& f, const double a, const double b);

/// HOOKE-JEEVES METHOD ///

constexpr double beta_init = 1.0;
constexpr double delta_init = 0.0;
constexpr double h_init = 0.01;
constexpr double h_lim = 1e-12;

void hooke_jeeves(const sigma_t& f, double& x0);
void hooke_jeeves_investigation(const sigma_t& f, const double x, double& h, double& x_out, bool reduce_step);
void hooke_jeeves(const sigma_t& f, double& x0, double& y0);
void hooke_jeeves_investigation(const sigma_t& f, const double x, const double y, double& hx, double& hy, double& x_out, double& y_out, bool reduce_step);

/// CONTRIBUTIONS ///

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

double sigma_kin(const double beta);
double sigma_cul(const double beta);
double sigma_x(const double beta);
double sigma_c(const double beta);
double sigma_kin2(const double beta);
double sigma_xc2(const double beta);
double sigma_kin4(const double beta);
double sigma_xc4(const double beta);
double sigma_ii(const double beta, const double delta);
double sigma_ei(const double beta, const double delta);

double D0(const double beta);
double D_ei(const double beta);
double D_delta(const double beta, const double delta);

const double mu = 0.5 * std::pow(3.0*M_PI*M_PI*n0, 2.0/3.0)
			      - std::pow(3.0*n0/M_PI, 1.0/3.0)
			      - (0.056*std::pow(n0, 2.0/3.0) + 0.0059*std::pow(n0, 1.0/3.0)) / 
			      	((0.079 + std::pow(n0, 1.0/3.0))*(0.079 + std::pow(n0, 1.0/3.0)))
			      - 0.4 * std::pow(Z, 2.0/3.0) * std::pow(4.0*M_PI*n0/3.0, 1.0/3.0)
			      + 4.0*M_PI*n0*r_c*r_c;



double sigma(const double beta, const double delta);
double work_function(const double beta, const double delta);

///

double f1(const double x, const double y)
{
	return x*x*x - x*x;
}


double f2(const double x, const double y)
{
	return (x-8)*(x-8) + (y-0.8)*(y-0.8);
}


int main()
{
	double x = 0.0, y = 0.0;
	hooke_jeeves(f2, x, y);
	printf("x = %5.6lf, y = %5.6lf\n", x, y);
}


double sigma(const double beta, const double delta)
{
	using std::sqrt, std::exp, std::sqrt;

	return sigma_kin(beta) + sigma_cul(beta) + sigma_x(beta) + sigma_c(beta)
		 + sigma_kin2(beta) + sigma_xc2(beta)
		 + sigma_kin4(beta) + sigma_xc4(beta)
		 + sigma_ii(beta, delta) + sigma_ei(beta, delta);
}


double work_function(const double beta, const double delta)
{
	return D0(beta) - mu + D_ei(beta) + D_delta(beta, delta);
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

	return (1.336 * pow(n0, 1.0/3.0) / (2160 * pow(3.0 * M_PI*M_PI, 2.0/3.0)))
	     * (
	     	pow(1 - 0.5*x, -11.0/3.0) * ((1 - 0.5*x)*(1 + x/16.0) + x*x/3.0)
	      + 29.0 * x*x * pow(0.5*x, -11.0/3.0) / 96.0
	     ) * x; 
}


double w_xc4(const double x)
{
	using std::pow, std::exp;

	return 7.35 * 1e-5 * (
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
	const w_t& f, 
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


double sigma_kin(const double beta)
{
	return w.kin / beta;
}

double sigma_cul(const double beta)
{
	return w.cul / (beta*beta*beta);
}

double sigma_x(const double beta)
{
	return w.x / beta;
}

double sigma_c(const double beta)
{
	return w.c / beta;
}

double sigma_kin2(const double beta)
{
	return w.kin2 * beta;
}

double sigma_xc2(const double beta)
{
	return w.xc2 * beta;
}

double sigma_kin4(const double beta)
{
	return w.kin4 * (beta*beta*beta);
}

double sigma_xc4(const double beta)
{
	return w.xc4 * (beta*beta*beta);
}

double sigma_ii(const double beta, const double delta)
{
	using std::sqrt, std::exp;
	const double coeff = sqrt(3.0)*Z*Z / (C*C*C);
	const double exp_arg = -4.0*M_PI*(d - 2.0*delta) / (sqrt(3.0)*C);
	
	return coeff * exp(exp_arg);
}

double sigma_ei(const double beta, const double delta)
{
	const double ei = 2.0*M_PI*n0*n0 * (
		       				1.0 - cosh(beta*r_c)*beta*d*exp(-0.5*beta*d) / (1.0 - exp(-beta*d))
		   			  ) / (beta*beta*beta);

	const double d_ei = 2.0 * M_PI * n0*n0 * d * (
		 					delta*delta + (1.0 - exp(beta*delta)) * exp(-0.5*beta*d) * cosh(beta*r_c) / (beta*beta)
		 				);

	return ei + d_ei; 
}


double D0(const double beta) 
{
	return 4.0*M_PI*n0 / (beta*beta);
}

double D_ei(const double beta)
{
	using std::exp, std::sinh, std::cosh;

	return -(4.0*M_PI*n0/(beta*beta)) * exp(-0.5*beta*d) * (
				beta*d*cosh(beta*r_c) - 2.0*sinh(-0.5*beta*d)
			) / (2.0 - exp(-beta*d));
}

double D_delta(const double beta, const double delta)
{
	using std::exp, std::sinh;

	return -4.0*M_PI*n0*exp(beta*(delta - 0.5*d))/(beta*beta*(2.0 - exp(-beta*d))) * (
				2.0 * (1.0 - exp(beta*delta)) * sinh(0.5*beta*d) - beta*beta*delta*d
			) + D_ei(beta) * (exp(beta*delta) - 1.0) - 4*M_PI*n0*d*delta;
}


void hooke_jeeves(const sigma_t& f, double& x0)
{
	double h = h_init;
	double x_prev = x0, x_next;
	bool reduce_step;

	while (h > h_lim) {
		reduce_step = true;
		hooke_jeeves_investigation(f, x_prev, h, x_next, reduce_step);

		double new_x_prev = x_prev + 2 * (x_next - x_prev);
		double new_x_next;

		reduce_step = false;
		hooke_jeeves_investigation(f, new_x_prev, h, new_x_next, reduce_step);

		while (f(new_x_next, 0.0) < f(x_next, 0.0)) {
			x_prev = x_next;
			x_next = new_x_next;

			new_x_prev = x_prev + 2 * (x_next - x_prev);
			hooke_jeeves_investigation(f, new_x_prev, h, new_x_next, reduce_step);
		}

		h /= 10;
	}

	x0 = x_next;
}


void hooke_jeeves_investigation(
	const sigma_t& f, const double x, 
	double& h, double& x_out, bool reduce_step
)
{
	x_out = x;

	if (f(x + h, 0.0) < f(x, 0.0)) {
		x_out = x + h;
	} else if (f(x - h, 0.0) < f(x, 0.0)) {
		x_out = x - h;
	} else if (reduce_step) {
		h /= 10;
	}
}


void hooke_jeeves(const sigma_t& f, double& x0, double& y0)
{
	double hx = h_init, hy = h_init;
	double x_prev = x0, y_prev = y0;
	double x_next, y_next;
	bool reduce_step;

	while (hx > h_lim && hy > h_lim) {
		reduce_step = true;
		hooke_jeeves_investigation(f, x_prev, y_prev, hx, hy, x_next, y_next, reduce_step);

		double new_x_prev = x_prev + 2.0 * (x_next - x_prev);
		double new_y_prev = y_prev + 2.0 * (y_next - y_prev);

		double new_x_next, new_y_next;

		reduce_step = false;
		hooke_jeeves_investigation(f, new_x_prev, new_y_prev, hx, hy, new_x_next, new_y_next, reduce_step);

		while (f(new_x_next, new_y_next) < f(x_next, y_next)) {
			x_prev = x_next;
			y_prev = y_next;

			x_next = new_x_next;
			y_next = new_y_next;

			new_x_prev = x_prev + 2.0 * (x_next - x_prev);
			new_y_prev = y_prev + 2.0 * (y_next - y_prev);

			hooke_jeeves_investigation(f, new_x_prev, new_y_prev, hx, hy, new_x_next, new_y_next, reduce_step);
		}
	}

	x0 = x_next;
	y0 = y_next;
}



void hooke_jeeves_investigation(
	const sigma_t& f, const double x, 
	const double y, double& hx,
	double& hy, double& x_out, 
	double& y_out, bool reduce_step
)
{
	x_out = x; y_out = y;

	if (f(x + hx, y) < f(x, y)) {
		x_out = x + hx;
	} else if (f(x - hx, y) < f(x, y)) {
		x_out = x - hx;
	} else if (reduce_step) {
		hx /= 10;
	} 

	if (f(x, y + hy) < f(x, y)) {
		y_out = y + hy;
	} else if (f(x, y - hy) < f(x, y)) {
		y_out = y - hy;
	} else if (reduce_step){
		hy /= 10;
	}
}