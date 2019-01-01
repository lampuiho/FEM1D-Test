#pragma once
#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/base/numbers.h>

class FEM1DGauSrc : public dealii::Function<1>
{
	// for filtering
	dealii::Point<1> location;
	double x1, y1, sigmax2, delay, sigma2;

public:
	FEM1DGauSrc(dealii::Point<1> location, double radius, double delay, double bandwidth);
	virtual void set_time(double new_time) final override;
	virtual double value(const dealii::Point<1> &p, const unsigned int  component = 0) const;
};
