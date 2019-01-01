#include "FEM1DGauSrc.h"

using namespace dealii;

FEM1DGauSrc::FEM1DGauSrc(Point<1> location, double radius, double delay, double bandwidth) : location(location), delay(delay)
{
	sigma2 = 4 * std::log(2) / std::pow(2 * numbers::PI *(bandwidth), 2);
	sigmax2 = std::pow(radius,2) / std::log(2);
}

void FEM1DGauSrc::set_time(double new_time)
{
	Function<1>::set_time(new_time);
	double x0 = exp(-std::pow((new_time - delay),2) / sigma2);
	double y0 = x0 - x1 + 0.995 * y1;
	x1 = x0;
	y1 = y0;
}

double FEM1DGauSrc::value(const Point<1>& p, const unsigned int component) const
{
	(void)component;
	Assert(component == 0, ExcIndexRange(component, 0, 1));

	double r2 = p.distance_square(this->location);
	return y1*exp(-r2 / sigmax2);
}
