#include <field/Field.h>

#include <field/Constant.h>





SHCoefficientFieldArray::SHCoefficientFieldArray( int order ):
	m_order(order)
{
	// initialize all coefficient fields to zero
	// NB: this is including coefficients with (l+m)%2==1
	// so the coefficient indices of PNSystem for 2d dont match those here,
	// which are the ones for general 3d
	int numCoeffs = (order + 1) * (order + 1);
	for( int i=0;i<numCoeffs;++i )
		m_coeff_fields.push_back(std::make_shared<Constant>(0.0));
}

void SHCoefficientFieldArray::setField( int l, int m, Field::Ptr field )
{
	// NB: this computes the sh index from l,m for the general 3d case
	m_coeff_fields[getSHIndex(l,m)] = field;
}

Field::Ptr SHCoefficientFieldArray::getField(int l, int m)
{
	// NB: this computes the sh index from l,m for the general 3d case
	return m_coeff_fields[getSHIndex(l,m)];
}

std::complex<double> SHCoefficientFieldArray::eval( int l, int m, const P3d& pWS )const
{
	// map l, m to shindex
	int shindex = l * (l + 1) + m; // NB: this computes the sh index from l,m for the general 3d case
	// evaluate field
	return m_coeff_fields[shindex]->eval(pWS);
}

int SHCoefficientFieldArray::getOrder() const
{
	return m_order;
}

int SHCoefficientFieldArray::getSHIndex(int l, int m) const
{
	return l * (l + 1) + m;
}

