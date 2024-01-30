// TableFunc.cpp
//
//////////////////////////////////////////////////////////////////////

#include "TableFunc.h"
#include <math.h>
#include "TableReader.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
CTableFunction::CTableFunction(CVector& _x, CVector& _y, double _leftVal, double _rightVal):
m_x(_x),
m_y(_y),
m_leftVal(_leftVal),
m_rightVal(_rightVal),
m_isReady(false),
m_ownsVectors(false),
m_table(0)
{

}

CTableFunction::CTableFunction(CTableReader* _table, double _leftVal, double _rightVal):
m_x(*(new CVector(_table->getColumn(0)))),
m_y(*(new CVector(_table->getColumn(1)))),
m_leftVal(_leftVal),
m_rightVal(_rightVal),
m_isReady(false),
m_ownsVectors(_table!=NULL),
m_table(_table)
{

}

CTableFunction::CTableFunction(const CRanges& _ranges, CVector& _y):
m_x(_ranges.m_midEnergies),
m_y(_y),
m_leftVal(0),
m_rightVal(0),
m_isReady(false),
m_ownsVectors(false),
m_table(0)
{

}

CTableFunction::~CTableFunction()
{
	delete m_table;
	if (m_ownsVectors) {
		delete &m_x;
		delete &m_y;
	}
}

void CTableFunction::init()
{
	if(!m_isReady)
	{
		prepare();
		m_isReady = true;
	}
}

void CTableFunction::takeDataOwnership(){
	m_ownsVectors=true;
}

namespace PropagDepreciated
{

void CLogScaleFunc::prepare()
{

	m_alphas.create(m_x.length());
	ASSERT(m_x.length()==m_y.length());
	int i;
	for(i=1;i<m_x.length();i++)
	{
		double y1=m_y[i-1];
		if(y1<=0.)
			continue;
		double x1=m_x[i-1];
		ASSERT(x1>0.);
		double y2=m_y[i];
		if(y2<=0.)
			continue;
		double x2=m_x[i];
		ASSERT(x2>0.);
		m_alphas[i-1] = log(y2/y1)/log(x2/x1);
	}
	m_alphas[i-1] = m_alphas[i-2];
}

CLogScaleFunc::~CLogScaleFunc()
{
 
}

double CLogScaleFunc::f(double _x) const
{
	ASSERT(_x>0.);
	if(m_x.length()<1)
		return 0.;
	int i = m_x.findLeftX(_x);
	if(i<0)
		return (_x < m_x[0])?m_leftVal:m_rightVal;
	double y1=m_y[i];
	double y2=m_y[i+1];
	if((y1<=0)&&(y2<=0.))
		return y1+(y2-y1)/(m_x[i+1]-m_x[i])*(_x - m_x[i]);
	return y1*pow(_x/m_x[i],m_alphas[i]);
}

/// also approximates left and right values
double CLogScaleFunc::fApprox(double _x) const
{
	ASSERT(_x>0.);
	if(m_x.length()<1)
		return 0.;
	int i = m_x.findLeftX(_x);
	if(i<0)
	{
		i = (_x < m_x[0])?0:(m_x.length()-2);
	};
	double y1=m_y[i];
	double y2=m_y[i+1];
	if((y1<=0)&&(y2<=0.))
		return y1+(y2-y1)/(m_x[i+1]-m_x[i])*(_x - m_x[i]);
	return y1*pow(_x/m_x[i],m_alphas[i]);
}

}//end PropagDepreciated

CCSplineFunc::~CCSplineFunc()
{
	if(m_spline) gsl_spline_free (m_spline);
	if(m_accel) gsl_interp_accel_free(m_accel);
}

double CCSplineFunc::f(double aX) const
{
	double result = (aX>=m_Xmin)?(aX<=m_Xmax?gsl_spline_eval(m_spline, aX, m_accel):m_rightVal):m_leftVal;
	return result;
	//y_deriv = gsl_spline_eval_deriv (m_spline, aX, m_accel);
	//y_deriv2 = gsl_spline_eval_deriv2 (m_spline, aX, m_accel);
}

double CCSplineFunc::fApprox(double aX) const
{
	double result = gsl_spline_eval(m_spline, aX, m_accel);
	return result;
}

void CCSplineFunc::prepare()
{
	m_accel = NULL;
	m_spline = NULL;
	m_accel = gsl_interp_accel_alloc ();
	m_spline = gsl_spline_alloc(gsl_interp_cspline, m_x.length());
	gsl_spline_init (m_spline, m_x.ptr(), m_y.ptr(), m_x.length());
	m_Xmin = m_x[0];
	m_Xmax = m_x[m_x.length()-1];
}

void CLogCSplineFunc::assertData()
{
	int length = m_x.length();
	int i=0;
	for(; i<length && m_y[i]==0.; i++);//skipping beginning zeros
	for(; i<length && m_y[i]>0.; i++);//skipping positive values
	for(; i<length && m_y[i]==0.; i++);//skipping ending zeros
	if (i<length) {
		ASSERT(false);
		cerr << "CLogCSplineFunc: incompatible data format" << endl;
	}
}

void CLogCSplineFunc::prepare()
{
	DEBUG_ONLY(assertData());
	m_accel = NULL;
	m_spline = NULL;
	int length = m_x.length();
	for(int i=0; i<length; i++){
		if(m_y[i]>0){//zero and negative point are thrown away!!!
			m_adjustedX.add(m_x[i]);
			m_logF.add(log(m_y[i]));
		}
	}
	m_accel = gsl_interp_accel_alloc ();
	m_spline = gsl_spline_alloc(gsl_interp_cspline, m_adjustedX.length());
	gsl_spline_init (m_spline, m_adjustedX.ptr(), m_logF.ptr(), m_adjustedX.length());
	m_Xmin = m_adjustedX[0];
	m_Xmax = m_adjustedX[m_adjustedX.length()-1];
}

double CLogCSplineFunc::f(double aX) const
{
	double result = (aX>=m_Xmin)?(aX<=m_Xmax?exp(gsl_spline_eval(m_spline, aX, m_accel)):m_rightVal):m_leftVal;
	return result;
	//y_deriv = gsl_spline_eval_deriv (m_spline, aX, m_accel);
	//y_deriv2 = gsl_spline_eval_deriv2 (m_spline, aX, m_accel);
}

double CLogCSplineFunc::fApprox(double aX) const
{
	double result = Exp(gsl_spline_eval(m_spline, aX, m_accel));
	return result;
}

CTableFunction* CCSplineFunc::getDerivative(int /*period*/) const
{
	int length = m_x.length();
	CVector* derX = new CVector(m_x);
	CVector* derY = new CVector(length);
	for(int i=0; i<length; i++){
		(*derY)[i] = gsl_spline_eval_deriv(m_spline, (*derX)[i], m_accel);
	}
	CTableFunction* result = createInstance(*derX, *derY);
	result->takeDataOwnership();
	return result;
}

CTableFunction* CCSplineFunc::get2ndDerivative() const
{
	int length = m_x.length();
	CVector* derX = new CVector(m_x);
	CVector* derY = new CVector(length);
	for(int i=0; i<length; i++){
		(*derY)[i] = gsl_spline_eval_deriv2(m_spline, (*derX)[i], m_accel);
	}
	CTableFunction* result = createInstance(*derX, *derY);
	result->takeDataOwnership();
	return result;
}

void CLinearFunc::prepare()
{
	ASSERT(m_x.length()>=2);
}

void CLinearFunc::SetExtension(EExtensionRule aLimits, bool aLeft, bool aRight)
{
	if(aLeft)
	{
		m_leftVal = 0;
		m_leftDerivative = 0.;
		if(aLimits != ExtNone)
		{
			m_leftVal = m_y[0];
			if(aLimits == ExtLinear)
				m_leftDerivative = (m_y[1]-m_leftVal)/(m_x[1]-m_x[0]);
		}

	}
	if(aRight)
	{
		m_rightVal = 0;
		m_rightDerivative = 0.;
		if(aLimits != ExtNone)
		{
			int last = m_y.length()-1;
			m_rightVal = m_y[last];
			if(aLimits == ExtLinear)
				m_rightDerivative = (m_rightVal-m_y[last-1])/(m_x[last]-m_x[last-1]);
		}
	}
}

double CLinearFunc::f(double _x) const
{
	int i = m_x.findLeftX(_x);
	
	if(i<0)
		return (_x < m_x[0]) ? (m_leftVal+(_x-m_x[0])*m_leftDerivative) : (m_rightVal+(_x-m_x[m_x.length()-1])*m_rightDerivative);

	double x1 = m_x[i];
	double y1 = m_y[i];
	double y2 = m_y[i+1];
	double x2 = m_x[i+1];

	return y1+(y2-y1)/(x2-x1)*(_x - x1);
}

double CLinearFunc::fApprox(double aX) const
{
	int i = m_x.findLeftX(aX);
	
	if(i<0){
		i=(aX < m_x[0])?0:m_x.length()-2;
	}

	double x1 = m_x[i];
	double y1 = m_y[i];
	double y2 = m_y[i+1];
	double x2 = m_x[i+1];

	return y1+(y2-y1)/(x2-x1)*(aX - x1);	
}

void CLogScaleLinearFunc::SetExtension(EExtensionRule aLimits, bool aLeft, bool aRight)
{
	if(aLimits==ExtNone){
		if(aLeft)
			fLogXY->setLeftVal(Log0);
		if(aRight)
			fLogXY->setRightVal(Log0);
	}
	else
		fLogXY->SetExtension(aLimits, aLeft, aRight);
}

double CLogScaleLinearFunc::f(double _x) const
{
	ASSERT(_x>0);
	return exp(fLogXY->f(log(_x) ));
}

double CLogScaleLinearFunc::fApprox(double _x) const
{
	ASSERT(_x>0);
	return exp(fLogXY->fApprox(log(_x) ));
}

const double CLogScaleLinearFunc::Log0 = -1000;//exp(-690)~1e-300 - close to minimal double

void CLogScaleLinearFunc::prepare()
{
	int iMax = m_x.length();
	ASSERT(iMax>=2);
	CVector* logX = new CVector(iMax);
	CVector* logY = new CVector(iMax);

	for(int i=0;i<iMax;i++)
	{
		ASSERT(m_x[i]>0);
		ASSERT(m_y[i]>0);
		(*logX)[i] = log(m_x[i]);
		(*logY)[i] = log(m_y[i]);
	}
	ASSERT(m_leftVal>=0 && m_rightVal>=0);
	double logLeftVal = m_leftVal>0 ? log(m_leftVal) : Log0;
	double logRightVal = m_rightVal>0 ? log(m_rightVal) : Log0;

	fLogXY = new CLinearFunc(*logX, *logY, logLeftVal, logRightVal);
	fLogXY->takeDataOwnership();
}

CTableFunction* CTableFunction::getDerivative(int period) const
{
	int iMax = (m_x.length()-1)/period;
	ASSERT(iMax>0);
	CVector* derX = new CVector(iMax);
	CVector* derY = new CVector(iMax);

	for(int i=0;i<iMax;i++)
	{	
		int i0 = i*period;
		int i1 = i0 + period;

		double dx = m_x[i1]-m_x[i0];
		ASSERT(dx>0);
		(*derX)[i] = 0.5*(m_x[i0]+m_x[i1]);
		(*derY)[i] = (m_y[i1]-m_y[i0])/dx;
	}
	CTableFunction* result = createInstance(*derX, *derY);
	result->takeDataOwnership();
	return result;
}

double CTableFunction::integrate(double g(double)) const
{
	CFunc f(g);
	return integrate(f);
}

double CTableFunction::integrate(const IFunction& g) const{
	/*
	int iMax = m_x.length()-1;
	double sum = 0;

	for(int i=0; i<iMax; i++)
	{
		double dx = m_x[i+1]-m_x[i];
		ASSERT(dx>0);
		double xMean = (m_x[i]+m_x[i+1])*0.5;
		double gMean = g(xMean);
		double yMean = 0.5*(m_y[i]+m_y[i+1]);
		
		sum += gMean*yMean*dx;
	}
	return sum;	*/

	int iMax = m_x.length();
	double sum = 0.;
	double product=m_y[0]*g(m_x[0]);
	for(int i=1;i<iMax;i++){
		double delta = product;
		product = m_y[i]*g(m_x[i]);
		delta += product;
		sum += delta*(m_x[i]-m_x[i-1]);
	}
	sum *= 0.5;
	return sum;
}

void CTableFunction::print(const char* fileName) const
{
	int iMax = m_x.length();
	CFilePtr file(Fopen(fileName, "wt"));
	for(int i=0; i<iMax; i++)
		fprintf(file, "%lg\t%lg\n", m_x[i], m_y[i]);
}

IFunction* CTableFunction::integral(double g(double)) const
{
	CFunc f(g);
	return integral(f);
}

IFunction* CTableFunction::integral(const IFunction& aFunc) const
{
	CVector* resultX = new CVector(m_x);
	CVector* resultY = new CVector(m_x.length());
	int iMax = m_x.length();
	double sum = 0.;
	(*resultY)[0] = 0.;
	double product=m_y[0]*aFunc(m_x[0]);
	for(int i=1;i<iMax;i++){
		double delta = product;
		product = m_y[i]*aFunc(m_x[i]);
		delta += product;
		sum += 0.5*delta*(m_x[i]-m_x[i-1]);
		(*resultY)[i] = sum;
	}
	CTableFunction* result = createInstance(*resultX, *resultY, 0, sum);
	result->takeDataOwnership();
	return result;
}

MatrixFunction::MatrixFunction(const std::string& aFile)
{
	std::ifstream dataFile(aFile.c_str());
	if(!dataFile)
		ThrowError("Failed to open " + aFile);
	std::string header;
	std::getline(dataFile, header);
	int logscX,logscY;
	const char* headerFormat = "# minX=%lg stepX=%lg logscaleX=%d minY=%lg stepY=%lg logscaleY=%d";
	if(sscanf(header.c_str(),headerFormat,
			&xMin,&xStep,&logscX,&yMin,&yStep,&logscY)!=6)
		ThrowError("Invalid header format in " + aFile + "\nExpected format: " + headerFormat);
	logScaleX = (bool)logscX;
	logScaleY = (bool)logscY;
	if((logScaleX && (xMin<=0 || xStep<=1.))||xStep<=0.)
		ThrowError("Invalid X scale in " + aFile);
	if((logScaleY && (yMin<=0 || yStep<=1.))||yStep<=0.)
			ThrowError("Invalid Y scale in " + aFile);
	fData = new CMatrix(dataFile);
	if(logScaleX) {
		xMax = xMin*pow(xStep,(*fData)[0].length()-1);
		xStep = log(xStep);//convert multiplier to log step
	}
	else{
		xMax = xMin + xStep*((*fData)[0].length()-1);
	}

	if(logScaleY){
		yMax = yMin*pow(yStep,fData->SizeI()-1);
		yStep = log(yStep);//convert multiplier to log step
	}
	else{
		yMax = yMin + yStep*(fData->SizeI()-1);
	}
}

double MatrixFunction::f(double x, double y) const
{
	double binX, binY;
	if(logScaleX)
	{//log scale
		binX = log(x/xMin)/xStep;
	}
	else
	{//linear scale
		binX = (x-xMin)/xStep;
	}
	if(binX<0. || binX>(*fData)[0].length()-1)
		return 0.;
	if(logScaleY)
	{//log scale
		binY = log(y/yMin)/yStep;
	}
	else
	{//linear scale
		binY = (y-yMin)/yStep;
	}
	if(binY<0. || binY>fData->SizeI()-1)
		return 0.;
	int iX = floor(binX);
	int iY = floor(binY);
	double weightX = 1.-binX+iX;
	double weightY = 1.-binY+iY;

	//linear f interpolation //TODO implement logscale f interpolation
	double f_11 = (*fData)[iY][iX]*weightX*weightY;
	double f_21 = (weightX<1.)?((*fData)[iY][iX+1]*(1.-weightX)*weightY) : 0.;
	double f_12 = (weightY<1.)?((*fData)[iY+1][iX]*weightX*(1.-weightY)) : 0.;
	double f_22 = (weightX<1.&&weightY<1.)?((*fData)[iY+1][iX+1]*(1.-weightX)*(1.-weightY)) : 0.;
	double result = f_11 + f_21 + f_12 + f_22;
	return result;
}
