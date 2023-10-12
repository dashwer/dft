#include <cstdio>
#include <memory>

#include <cmath>

constexpr int accuracy = 3; 											  
constexpr double order = std::pow(10, accuracy);


constexpr double h = std::round(6.62607015 * order) / order * 1e-34;       ///Постоянная Планка("6.62607015e-34")
constexpr double c = std::round(2.99792458 * order) / order * 1e8;         ///Скорость света("299792458")
constexpr double e = std::round(1.602176634 * order) / order * 1e-19;      ///Заряд электрона("1.602176634e-19")
constexpr double m_e = std::round(9.1093837015 * order) / order * 1e-31;   ///Масса электрона("9.1093837015e-31")
constexpr double pi = std::round(3.14159265358979 * order) / order; 	   ///Число пи("3.14159265358979")


constexpr double a_0 = (h*h * 1e7) / (m_e * 4 * pi*pi * e*e * c*c);        /// Боровский радиус
constexpr double E_h = (m_e * 4*pi*pi * e*e*e*e * c*c*c*c * 1e-14) / (h*h);/// Энергия Хартри 

/*
measures_t - список единиц измерений
*/
enum class measures_t {
	Bohrs, Meters, Angstroms, 
	Joules, eV, Hartree,
	Joules_meters, Hartree_Bohrs,
	Newton, Atomic_force_unit
};

/*--------------------------------------------------------
Базовый класс физической величины:

Поля:
value - значение физической величины
measure - единица измерения

Методы:
get_value - получить значение физической величины 
set_value - установить значений физической величины
get_measure - получить единицу измерения
set_measure - установить единицу измерения
printer - форматированный вывод значения физической 
		  величины
converter - конвертирует значение исходя из 
            единиц измерения
get_value_impl - реализация метода get_value
set_value_impl - реализация метода set_value
get_measure_impl - реализация метода get_measure
set_measure_impl - реализация метода set_measure
printer_impl - реализация метода printer
--------------------------------------------------------*/
class Physical_quantity_t
{
public:
	Physical_quantity_t() = default;

	virtual ~Physical_quantity_t() = default;

	double get_value() const
	{
		return this->get_value_impl();
	}

	int set_value(double input_value, measures_t input_measure)
	{
		return this->set_value_impl(input_value, input_measure);
	}

	measures_t get_measure() const
	{
		return this->get_measure_impl();
	}

	int set_measure(measures_t input_measure)
	{
		return this->set_measure_impl(input_measure);
	}

	int printer() const
	{
		return this->printer_impl();
	}

protected:
	double value;
	measures_t measure;

	virtual double converter(double input_value, measures_t input_measure, measures_t use_measure) const = 0;
	virtual double get_value_impl() const = 0;
	virtual int set_value_impl(double input_value, measures_t input_measure) = 0;
	virtual measures_t get_measure_impl() const = 0;
	virtual int set_measure_impl(measures_t input_measure) = 0;
	virtual int printer_impl() const = 0;
};

/*--------------------------------------------------------
Производный класс единиц измерения длины;
--------------------------------------------------------*/
class Distance_t : public Physical_quantity_t
{
public:
	Distance_t(double input_value, measures_t input_measure, measures_t use_measure)
	{
		this->measure = use_measure;
		this->value = this->converter(input_value, use_measure, input_measure);
	}

	Distance_t(std::unique_ptr<Distance_t> &obj, measures_t use_measure)
	{
		this->measure = use_measure;
		this->value = this->converter(obj->get_value(), use_measure, obj->get_measure());
	}	

private:
	virtual double converter(double input_value, measures_t input_measure, measures_t use_measure) const override
	{
		switch(use_measure)
		{
			case measures_t::Bohrs:
				switch(input_measure)
				{
					case measures_t::Meters:
						return input_value * a_0;
						break;
					case measures_t::Angstroms:
						return input_value * a_0 * 1e10;
						break;
					case measures_t::Bohrs:
						return input_value;
						break;
					default:
						return -1.0;
				}
				break;

			case measures_t::Angstroms:
				switch(input_measure)
				{
					case measures_t::Meters:
						return input_value * 1e-10;
						break;
					case measures_t::Angstroms:
						return input_value;
						break;
					case measures_t::Bohrs:
						return input_value * 1e-10 / a_0;
						break;
					default:
						return -1.0;
				}
				break;

			case measures_t::Meters:
				switch(input_measure)
				{
					case measures_t::Meters:
						return input_value;
						break;
					case measures_t::Angstroms:
						return input_value * 1e10;
						break;
					case measures_t::Bohrs:
						return input_value / a_0;
						break;
					default:
						return -1.0;
				}
				break;

			default:
				return -1.0;
		}
	}

	virtual double get_value_impl() const override
	{
		return this->value;
	}

	virtual int set_value_impl(double input_value, measures_t input_measure) override
	{
		this->value = this->converter(input_value, this->measure, input_measure);
		return 0;
	}

	virtual measures_t get_measure_impl() const override
	{
		return this->measure;
	}

	virtual int set_measure_impl(measures_t input_measure) override
	{
		this->value = this->converter(this->value, input_measure, this->measure);
		this->measure = input_measure;
		return 0;
	}
	
	virtual int printer_impl() const override
	{
		switch(this->measure)
		{
			case measures_t::Bohrs:
				std::printf("Distance value = %.12e Bohrs \n", this->value);
				break;
			case measures_t::Meters:
				std::printf("Distance value = %.12e Meters \n", this->value);
				break;
			case measures_t::Angstroms:
				std::printf("Distance value = %.12e Angstroms \n", this->value);
				break;
			default:
				std::printf("Incorrect call \n");
		}

		return 0;
	}
};

/*--------------------------------------------------------
Производный класс единиц измерения плотности;
--------------------------------------------------------*/
class Density_t : public Physical_quantity_t
{
public:
	Density_t(double input_value, measures_t input_measure, measures_t use_measure)
	{
		this->measure = use_measure;
		this->value = this->converter(input_value, use_measure, input_measure);
	}

	Density_t(std::unique_ptr<Physical_quantity_t> &obj, measures_t use_measure)
	{
		this->measure = use_measure;
		this->value = this->converter(obj->get_value(), use_measure, obj->get_measure());
	}	

private:
	virtual double converter(double input_value, measures_t input_measure, measures_t use_measure) const override
	{
		switch (use_measure)
		{
			case measures_t::Bohrs:

				switch (input_measure)
				{
					case measures_t::Meters:
						return input_value / (a_0 * a_0 * a_0);
						break;

					case measures_t::Bohrs:
						return input_value;
						break;

					default:
						return -1.0;			
				}
				break;
		
			case measures_t::Meters:

				switch (input_measure)
				{
					case measures_t::Meters:
						return input_value;
						break;

					case measures_t::Bohrs:
						return input_value * (a_0 * a_0 * a_0);
						break;

					default:
						return -1.0;			
				}
				break;
			
			default:
				return -1.0;
		}
		
	}

	virtual double get_value_impl() const override
	{
		return this->value;
	}

	virtual int set_value_impl(double input_value, measures_t input_measure) override
	{
		this->value = this->converter(input_value, this->measure, input_measure);
		return 0;
	}

	virtual measures_t get_measure_impl() const override
	{
		return this->measure;
	}

	virtual int set_measure_impl(measures_t input_measure) override
	{
		this->value = this->converter(this->value, input_measure, this->measure);
		this->measure = input_measure;
		return 0;
	}
	
	virtual int printer_impl() const override
	{
		if(this->measure == measures_t::Bohrs)
		{
			std::printf("Density value = %.12e Bohrs^-3 \n", this->value);
		}
		else if(this->measure == measures_t::Meters)
		{
			std::printf("Density value = %.12e Meters^-3 \n", this->value);
		}
		else
		{
			std::printf("Incorrect call \n");
		}
		return 0;
	}
};

/*--------------------------------------------------------
Производный класс единиц измерения энергии;
--------------------------------------------------------*/
class Energy_t : public Physical_quantity_t
{
public:
	Energy_t(double input_value, measures_t input_measure, measures_t use_measure)
	{
		this->measure = use_measure;
		this->value = this->converter(input_value, use_measure, input_measure);
	}

	Energy_t(std::unique_ptr<Physical_quantity_t> &obj, measures_t use_measure)
	{
		this->measure = use_measure;
		this->value = this->converter(obj->get_value(), use_measure, obj->get_measure());
	}

private:
	virtual double converter(double input_value, measures_t input_measure, measures_t use_measure) const override
	{
		switch(use_measure)
		{
			case measures_t::Joules:
				switch(input_measure)
				{
					case measures_t::Joules:
						return input_value;
						break;
					case measures_t::eV:
						return input_value / e;
						break;
					case measures_t::Hartree:
						return input_value / E_h;
						break;
					default:
						return -1.0;
				}
				break;

			case measures_t::eV:
				switch(input_measure)
				{
					case measures_t::Joules:
						return input_value * e;
						break;
					case measures_t::eV:
						return input_value;
						break;
					case measures_t::Hartree:
						return input_value * e / E_h;
						break;
					default:
						return -1.0;
				}
				break;

			case measures_t::Hartree:
				switch(input_measure)
				{
					case measures_t::Joules:
						return input_value * E_h;
						break;
					case measures_t::eV:
						return input_value * E_h / e;
						break;
					case measures_t::Hartree:
						return input_value;
						break;
					default:
						return -1.0;
				}
				break;

			default:
				return -1.0;
		}
	}

	virtual double get_value_impl() const override
	{
		return this->value;
	}

	virtual int set_value_impl(double input_value, measures_t input_measure) override
	{
		this->value = this->converter(input_value, this->measure, input_measure);
		return 0;
	}

	virtual measures_t get_measure_impl() const override
	{
		return this->measure;
	}

	virtual int set_measure_impl(measures_t input_measure) override
	{
		this->value = this->converter(this->value, input_measure, this->measure);
		this->measure = input_measure;
		return 0;
	}

	virtual int printer_impl() const override
	{
		switch(this->measure)
		{
			case measures_t::Joules:
				std::printf("Energy value = %.12e Joules \n", this->value);
				break;
			case measures_t::eV:
				std::printf("Energy value = %.12e eV \n", this->value);
				break;
			case measures_t::Hartree:
				std::printf("Energy value = %.12e Hartree \n", this->value);
				break;
			default:
				std::printf("Incorrect call \n");
		}
		return 0;
	}
};

/*--------------------------------------------------------
Производный класс единиц измерения поверхностной энергии;
--------------------------------------------------------*/
class Surface_energy_t : public Physical_quantity_t
{
public:
	Surface_energy_t(double input_value, measures_t input_measure, measures_t use_measure)
	{
		this->measure = use_measure;
		this->value = this->converter(input_value, use_measure, input_measure);
	}

	Surface_energy_t(std::unique_ptr<Physical_quantity_t> &obj, measures_t use_measure)
	{
		this->measure = use_measure;
		this->value = this->converter(obj->get_value(), use_measure, obj->get_measure());
	}

private:
	virtual double converter(double input_value, measures_t input_measure, measures_t use_measure) const override
	{
		switch(use_measure)
		{
			case measures_t::Joules_meters:
				switch(input_measure)
				{
					case measures_t::Joules_meters:
						return input_value;
						break;
					case measures_t::Hartree_Bohrs:
						return input_value / (E_h / (a_0 * a_0));
						break;
					default:
						return -1.0;
				}
				break;

			case measures_t::Hartree_Bohrs:
				switch(input_measure)
				{
					case measures_t::Joules_meters:
						return input_value * (E_h / (a_0 * a_0));
						break;
					case measures_t::Hartree_Bohrs:
						return input_value;
						break;
					default:
						return -1.0;
				}
				break;

			default:
				return -1.0;
		}
	}

	virtual double get_value_impl() const override
	{
		return this->value;
	}

	virtual int set_value_impl(double input_value, measures_t input_measure) override
	{
		this->value = this->converter(input_value, this->measure, input_measure);
		return 0;
	}

	virtual measures_t get_measure_impl() const override
	{
		return this->measure;
	}

	virtual int set_measure_impl(measures_t input_measure) override
	{
		this->value = this->converter(this->value, input_measure, this->measure);
		this->measure = input_measure;
		return 0;
	}

	virtual int printer_impl() const override
	{
		switch(this->measure)
		{
			case measures_t::Joules_meters:
				std::printf("Surface energy value = %.12e Joules/meters^2 \n", this->value);
				break;
			case measures_t::Hartree_Bohrs:
				std::printf("Surface energy value = %.12e Hartree/Bohrs^2 \n", this->value);
				break;
			default:
				std::printf("Incorrect call \n");
		}
		return 0;
	}
};

/*--------------------------------------------------------
Производный класс единиц измерения силы;
--------------------------------------------------------*/
class Force_t : public Physical_quantity_t
{
public:
	Force_t(double input_value, measures_t input_measure, measures_t use_measure)
	{
		this->measure = use_measure;
		this->value = this->converter(input_value, use_measure, input_measure);
	}

	Force_t(std::unique_ptr<Physical_quantity_t> &obj, measures_t use_measure)
	{
		this->measure = use_measure;
		this->value = this->converter(obj->get_value(), use_measure, obj->get_measure());
	}

private:
	virtual double converter(double input_value, measures_t input_measure, measures_t use_measure) const override
	{
		switch(use_measure)
		{
			case measures_t::Newton:
				switch(input_measure)
				{
					case measures_t::Newton:
						return input_value;
						break;
					case measures_t::Atomic_force_unit:
						return input_value / (E_h / a_0);
						break;
					default:
						return -1.0;
				}
				break;

			case measures_t::Atomic_force_unit:
				switch(input_measure)
				{
					case measures_t::Newton:
						return input_value * (E_h / a_0);
						break;
					case measures_t::Atomic_force_unit:
						return input_value;
						break;
					default:
						return -1.0;
				}
				break;

			default:
				return -1.0;
		}
	}

	virtual double get_value_impl() const override
	{
		return this->value;
	}

	virtual int set_value_impl(double input_value, measures_t input_measure) override
	{
		this->value = this->converter(input_value, this->measure, input_measure);
		return 0;
	}

	virtual measures_t get_measure_impl() const override
	{
		return this->measure;
	}

	virtual int set_measure_impl(measures_t input_measure) override
	{
		this->value = this->converter(this->value, input_measure, this->measure);
		this->measure = input_measure;
		return 0;
	}

	virtual int printer_impl() const override
	{
		switch(this->measure)
		{
			case measures_t::Newton:
				std::printf("Force value = %.12e Newtons \n", this->value);
				break;
			case measures_t::Atomic_force_unit:
				std::printf("Force value = %.12e a.u. \n", this->value);
				break;
			default:
				std::printf("Incorrect call \n");
		}
		return 0;
	}
};

/*--------------------------------------------------------
Функции для теста основных методов классов:
Distance_t, Density_t, Energy_t, Surface_energy_t, Force_t;
--------------------------------------------------------*/
void test_Dist_class()
{
	double input_value = 1.5e+22;
	std::unique_ptr<Distance_t> dist_3 = std::make_unique<Distance_t>(input_value, measures_t::Meters, measures_t::Angstroms);
	std::unique_ptr<Physical_quantity_t> dist_1 = std::make_unique<Distance_t>(input_value, measures_t::Meters, measures_t::Angstroms);
	std::unique_ptr<Physical_quantity_t> dist_2 = std::make_unique<Distance_t>(dist_3, measures_t::Bohrs);
	dist_1->printer();
	dist_1->set_measure(measures_t::Bohrs);
	dist_1->printer();

	dist_2->printer();
	input_value = 4.5e+12;
	dist_2->set_value(input_value, measures_t::Meters);
	dist_2->printer();
}

void test_Density_class()
{
	double input_value = 1.5e+12;
	std::unique_ptr<Physical_quantity_t> den_1 = std::make_unique<Density_t>(input_value, measures_t::Bohrs, measures_t::Meters);
	std::unique_ptr<Physical_quantity_t> den_2 = std::make_unique<Density_t>(den_1, measures_t::Bohrs);
	den_1->printer();
	den_1->set_measure(measures_t::Meters);
	den_1->printer();

	den_2->printer();
	input_value = 4.5e+12;
	den_2->set_value(input_value, measures_t::Meters);
	den_2->printer();
}

void test_Energy_class()
{
	double input_value = 1.5e+3;
	std::unique_ptr<Physical_quantity_t> E_1 = std::make_unique<Energy_t>(input_value, measures_t::Joules, measures_t::eV);
	std::unique_ptr<Physical_quantity_t> E_2 = std::make_unique<Energy_t>(E_1, measures_t::Hartree);
	E_1->printer();
	E_1->set_measure(measures_t::Hartree);
	E_1->printer();

	E_2->printer();
	input_value = 4.5e+12;
	E_2->set_value(input_value, measures_t::Joules);
	E_2->printer();
}

void test_Surface_energy_class()
{
	double input_value = 1.5e+12;
	std::unique_ptr<Physical_quantity_t> surf_1 = std::make_unique<Surface_energy_t>(input_value, measures_t::Hartree_Bohrs, measures_t::Joules_meters);
	std::unique_ptr<Physical_quantity_t> surf_2 = std::make_unique<Surface_energy_t>(surf_1, measures_t::Hartree_Bohrs);
	surf_1->printer();
	surf_1->set_measure(measures_t::Hartree_Bohrs);
	surf_1->printer();

	surf_2->printer();
	input_value = 4.5e+12;
	surf_2->set_value(input_value, measures_t::Joules_meters);
	surf_2->printer();
}


void test_Force_class()
{
	double input_value = 1.5e+12;
	std::unique_ptr<Physical_quantity_t> F_1 = std::make_unique<Force_t>(input_value, measures_t::Atomic_force_unit, measures_t::Newton);
	std::unique_ptr<Physical_quantity_t> F_2 = std::make_unique<Force_t>(F_1, measures_t::Atomic_force_unit);
	F_1->printer();
	F_1->set_measure(measures_t::Atomic_force_unit);
	F_1->printer();

	F_2->printer();
	input_value = 4.5e+12;
	F_2->set_value(input_value, measures_t::Newton);
	F_2->printer();
}


int main()
{
	std::unique_ptr<Surface_energy_t> dist = std::make_unique<Surface_energy_t>(
											1.0, measures_t::Hartree_Bohrs, 
											measures_t::Joules_meters
										);
	dist->printer();
}