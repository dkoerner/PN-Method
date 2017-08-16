#pragma once
#include <memory>
#include <common/Domain.h>
#include <common/GridLocation.h>



struct Field
{
	typedef std::shared_ptr<Field> Ptr;
	virtual std::complex<double> eval( const P2d& pWS )const=0;
	virtual std::complex<double> dx(const P2d& pWS)const=0;
	virtual std::complex<double> dxdx(const P2d& pWS)const=0;
	virtual std::complex<double> dxdy(const P2d& pWS)const=0;
	virtual std::complex<double> dy(const P2d& pWS)const=0;
	virtual std::complex<double> dydy(const P2d& pWS)const=0;
	virtual std::complex<double> dydx(const P2d& pWS)const=0;
	virtual std::complex<double> dz(const P2d& pWS)const=0;
};


struct RadianceField
{
	typedef std::shared_ptr<RadianceField> Ptr;
	virtual std::complex<double> eval( const P2d& pWS, const V3d& omega )const=0;
	virtual std::complex<double> dx(const P2d& pWS, const V3d& omega)const=0;
	virtual std::complex<double> dxdx(const P2d& pWS, const V3d& omega)const=0;
	virtual std::complex<double> dxdy(const P2d& pWS, const V3d& omega)const=0;
	virtual std::complex<double> dy(const P2d& pWS, const V3d& omega)const=0;
	virtual std::complex<double> dydy(const P2d& pWS, const V3d& omega)const=0;
	virtual std::complex<double> dydx(const P2d& pWS, const V3d& omega)const=0;
	virtual std::complex<double> dz(const P2d& pWS, const V3d& omega)const=0;
	virtual std::complex<double> integral_over_solid_angle(const P2d& pWS)const=0;
};










