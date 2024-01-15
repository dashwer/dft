#include <cstdio>
#include <cstdlib>
#include <cmath>

#define PHI_COMPUTE_ERROR -1.0

constexpr double minus_inf = -10.0;
constexpr double plus_inf = 15.0;

constexpr double D = 1.0;
constexpr double h = 5.0;
constexpr double eps = 1.0;
constexpr double n_0 = 0.02;
constexpr double beta = std::sqrt((6*M_PI*n_0) / (0.5*std::pow(3*M_PI*M_PI*n_0, 2.0/3.0)));
constexpr double k = -4.0 * M_PI * n_0 / (beta * beta); 

constexpr double delta = std::sqrt(eps);

constexpr double step = 0.01;
constexpr int steps_amount = (plus_inf - minus_inf) / step;

double C_2;
double C_3;
double C_4;
double C_5;
double C_6;
double C_7;

int const_config()
{
    C_3 = 0.5 * k;
    
    C_4 = ( 
        C_3 * std::exp(-2*beta*D) * (delta*delta - 1) 
            * (std::exp(beta*h/delta) - std::exp(-beta*h/delta))
        ) / (
          (delta*delta - 1)*(delta*delta - 1)*std::exp(-beta*h/delta)
          - (delta*delta + 1)*(delta*delta + 1)*std::exp(beta*h/delta)
    );

    C_5 = std::exp(beta*D/delta) * (
            (delta + 1) * C_3 * std::exp(-beta*D)
          + (delta - 1) * C_4 * std::exp(beta*D)
    ) / (2*delta);

    C_6 = std::exp(-beta*D/delta) * (
            (delta - 1) * C_3 * std::exp(-beta*D)
          + (delta + 1) * C_4 * std::exp(beta*D)
    ) / (2*delta);

    C_7 = std::exp(beta*(D+h)) * (
          C_5 * std::exp(-beta*(D+h)/delta)
        + C_6 * std::exp(beta*(D+h)/delta)
    ); 
    
    C_2 = C_4 - C_3;

    return 0;
}

double phi_1(double x)
{
    return C_2 * std::exp(beta * x) + k;
}

double phi_2(double x)
{
    return C_3 * std::exp(-beta * x) + C_4 * std::exp(beta * x);
}

double phi_3(double x)
{
    return C_5 * std::exp(-beta * x / delta) + C_6 * std::exp(beta * x / delta);
}

double phi_4(double x)
{
    return C_7 * std::exp(-beta * x);
}


double phi(double x)
{
    if(x <= 0.0)
    {
        return phi_1(x);
    }
    else if(x >= 0.0 && x <= D)
    {
        return phi_2(x);
    }
    else if(x >= D && x <= D + h)
    {
        return phi_3(x);
    }
    else if(x >= D + h)
    {
        return phi_4(x);
    }
    else
    {
        return PHI_COMPUTE_ERROR;
    }
}

double n(double x)
{
    return -beta * beta * phi(x) / (4 * M_PI);
}

int main()
{
    const_config();

    char file_name_phi[256];
    char file_name_n[256];

    std::sprintf(
        file_name_phi, 
        "data/phi[eps=%.3lf][beta=%.3lf][D=%.3lf].dat", 
        eps, beta, D
    );

    std::sprintf(
        file_name_n, 
        "data/n[eps=%.3lf][beta=%.3lf][D=%.3lf].dat", 
        eps, beta, D
    );

    FILE* out_full_phi = fopen(file_name_phi, "w");
    if(out_full_phi == nullptr)
    {
        exit(-1);
    }

    FILE* out_full_n = fopen(file_name_n, "w");
    if(out_full_n == nullptr)
    {
        exit(-1);
    }

    for(int step_num = 0; step_num <= steps_amount; step_num++)
    {
        std::fprintf(
            out_full_phi, "%lf\t%lf\n", 
            minus_inf + step * step_num, 
            phi(minus_inf + step * step_num)
        );

        std::fprintf(
            out_full_n, "%lf \t %lf \n", 
            minus_inf + step * step_num, 
            n(minus_inf + step * step_num)
        );
    }

    fclose(out_full_phi);
    fclose(out_full_n);
}