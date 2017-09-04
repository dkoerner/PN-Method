#pragma once
#include <memory>
#include <common/Domain.h>


// Abstract interface to a (complex-valued) scalarfield which
// can be evaluated at worldspace positions
struct Field
{
	typedef std::shared_ptr<Field> Ptr;
	virtual std::complex<double> eval( const P2d& pWS )const=0;
	/*
	virtual std::complex<double> dx(const P2d& pWS)const=0;
	virtual std::complex<double> dxdx(const P2d& pWS)const=0;
	virtual std::complex<double> dxdy(const P2d& pWS)const=0;
	virtual std::complex<double> dy(const P2d& pWS)const=0;
	virtual std::complex<double> dydy(const P2d& pWS)const=0;
	virtual std::complex<double> dydx(const P2d& pWS)const=0;
	virtual std::complex<double> dz(const P2d& pWS)const=0;
	*/
	virtual void test()const{}
};





// Some RTE parameters also depend on the SH band given by l and m. Therefore these
// parameters have a field per coefficient and constitute an array of fields. This is
// what this class is for.
struct SHCoefficientFieldArray
{
	typedef std::shared_ptr<SHCoefficientFieldArray> Ptr;



	SHCoefficientFieldArray( int order );

	void setField( int l, int m, Field::Ptr field );
	Field::Ptr getField(int l, int m);
	std::complex<double> eval( int l, int m, const V2d& pWS )const;

	int getOrder()const;

private:
	// This function computes the linear SH index from l,m indices
	// NB: this uses the indexing for the general 3d case
	int getSHIndex( int l, int m )const;

	std::vector<Field::Ptr> m_coeff_fields; // a field for each coefficient
	int m_order;
};



// Interface to a radiancefield, which depends on position and a direction.
// It is used for using the solution of a pn problem in applications
// and for validating PN terms for a given radiance field.
struct RadianceField
{
	typedef std::shared_ptr<RadianceField> Ptr;

	virtual std::complex<double> eval( const P2d& pWS, const V3d& omega )const=0;
	/*
	virtual std::complex<double> dx(const P2d& pWS, const V3d& omega)const=0;
	virtual std::complex<double> dxdx(const P2d& pWS, const V3d& omega)const=0;
	virtual std::complex<double> dxdy(const P2d& pWS, const V3d& omega)const=0;
	virtual std::complex<double> dy(const P2d& pWS, const V3d& omega)const=0;
	virtual std::complex<double> dydy(const P2d& pWS, const V3d& omega)const=0;
	virtual std::complex<double> dydx(const P2d& pWS, const V3d& omega)const=0;
	virtual std::complex<double> dz(const P2d& pWS, const V3d& omega)const=0;
	virtual std::complex<double> integral_over_solid_angle(const P2d& pWS)const=0;
	*/
};





