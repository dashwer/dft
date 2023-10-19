#include <array>
#include <cmath>
#include <functional>
#include <stdio.h>

constexpr int nodes_number = 8;
constexpr double beta_init = 0.8;
constexpr double h_init = 0.01;
constexpr double h_lim = 1e-8;
constexpr double n0 = 0.0194;

constexpr double weights[nodes_number] = 
{
	0.10122854,0.22238104,
	0.31370664,0.36268378,
	0.36268378,0.31370664,
	0.22238104,0.10122854
};
	     
constexpr double nodes[nodes_number] = 
{
	-0.96028986,-0.79666648,
	-0.52553242,-0.18343464,
	0.18343464,0.52553242,
	0.79666648,0.96028986
};

using func_t = std::function<double(const double)>;
using nodes_container_t = std::array<double, nodes_number>;

double w_kin(const double x);
double w_cul(const double x);
double w_x(const double x);
double w_c(const double x);
//double w_kin2(const double x);
double w_xc2(const double x);

struct
{
	double kin, cul, x, c;
	double kin2, xc2;
} w;

nodes_container_t nodes_recalculation(const double a, const double b);
double gauss_integral(const func_t& f);
double sigma(const double beta);
double hooke_jeeves(const func_t& f, const double x0);
double hooke_jeeves_investigation(const func_t& f, const double x, const double h);


int main()
{
	w.kin = gauss_integral(w_kin);
	w.cul = gauss_integral(w_cul);
	w.x = gauss_integral(w_x);
	w.c = gauss_integral(w_c);
	w.kin2 = n0 * std::log(2.0) / 72;
	w.xc2 = gauss_integral(w_xc2);

	const double beta0 = hooke_jeeves(sigma, beta_init);
	const double sigma0 = sigma(beta0);

	printf("Result: %5.8f\n", sigma0);
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
			   	 (1.0 - 0.5*x)*(1 - 0.5*x) + 0.25*x*x - 1.0
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
	return -0.0056 * std::pow(n0, 4.0/3.0) * (
				
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


double gauss_integral(const func_t& f)
{
	nodes_container_t recalculated_nodes = nodes_recalculation(0.0, 1.0);
	double I = 0.0;

	for (int i = 0; i < nodes_number; ++i) {
		I += weights[i] * f(recalculated_nodes[i]);
	}
	I /= 2.0;

	return I;
}


double sigma(const double beta)
{
	return w.kin / beta 
		 + w.cul / (beta*beta*beta)
		 + w.x / beta + w.c / beta
		 + w.kin2 * beta + w.xc2 * beta;
}


double hooke_jeeves(
	const func_t& f, 
	const double x0
)
{
	double h = h_init;
	double x_prev = x0, x_next;

	while (h > h_lim) {
		x_next = hooke_jeeves_investigation(f, x_prev, h);

		if (f(x_next) >= f(x_prev)) {	
			h /= 10;
			continue;
		}

		double new_x_prev = x_prev + 2 * (x_next - x_prev);
		double new_x_next = hooke_jeeves_investigation(f, new_x_prev, h);
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


double hooke_jeeves_investigation(
	const func_t& f, 
	const double x,
	const double h
)
{
	double new_x;

	if (f(x + h) < f(x)) {
		new_x = x + h;
	} else if (f(x - h) < f(x)) {
		new_x = x - h;
	} else {
		new_x = x;
	}

	return new_x;
}