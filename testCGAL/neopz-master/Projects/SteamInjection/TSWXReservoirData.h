//---------------------------------------------------------------------------

#ifndef TSWXReservoirDataH
#define TSWXReservoirDataH
//---------------------------------------------------------------------------

#include "ThermalMethodsTables.h"

class TSwxReservoirData
{

	public:

	TSwxReservoirData();

	void SetData(double rw, double re, double Hr, double Pe,
					double Kh, double Poros, double So, double Tr,
					double RockDensity, double RockSpecificHeat,
					double OilDensity, double OilSpecificHeat, double OilViscosity);

	// M�todos para determinar as propriedades do reservatorio, considerando a rocha e as saturacoes dos fluidos em alguma posicao radial

	// M�todo que retorna a densidade e calor espec�fico do reservat�rio (rocha e flu�dos) na regi�o aquecida
	// As satura��es da �gua, vapor e �leo � dado como m�dia na regi�o aquecida (de vapor)
	// A satura��o do vapor � dependente da satura��o da �gua e �leo: Sv = 1 - (Sw + So)
	double Rho_C_Estrela(WaterDataInStateOfSaturation &table, double temp, double Sw, double So);   // Unidades: [kg/m3][kJ/(kg K)] --> [kJ/(m3 K)]

//	void setReservoirTemperature(double temp) { fTo = temp; }   // Unidades: C (graus Celsius),  K (graus Kelvin)
//
//	// M�todos para recuperar dados do reservat�rio
	double Temperature() { return fTr; }
	double getHReservoir() { return fHr; }
	double getDrainageRadius() {return fre; }
	double getWellRadius() { return frw; }
	double getStaticPressure() { return fPe; }
	double getKh() {return fKh; }
	double OilSaturation() { return fSo; }
	double OilDensity() { return fOilDensity; }
	double OilVisc() { return fOilVisc; }
	double OilSpecificHeat() { return fOilSpecificHeat; }
	double getPorosity() { return fPoros; }
	double getRockDensity() { return fRockDensity; }
	double getRockSpecificHeat() { return fRockSpecificHeat; }

	protected:
	std::vector< std::vector<double> > fTable;

	/**
	 * frw: raio do poco [m]
	 * fre: raio de drenagem [m]
	 * fHr: espessura do reservatorio [m]
	 * fPe: pressao estatica [Pa]
	 * fKh: permeabilidade horizontal [m2]
	 * fPoros: porosidade (0~1)
	 * fSo: saturacao de Oleo (0~1)
	 * fTr: temperatura inicial (grau Celsius)
	 * fRockDensity: Massa especifica [kg/m3]
	 * fRockSpecificHeat: Calor especifico [kJ/(kg C)]
	 *
	 * Oleo:
	 * fOilDensity: Massa especifica [kg/m3]
	 * fOilSpecificHeat: Calor especifico [kJ/(kg C)]
	 * fOilVisc: Viscosidade
	 */
	double frw, fre, fHr, fPe, fKh, fPoros, fSo, fTr, fRockDensity, fRockSpecificHeat, fOilDensity, fOilSpecificHeat, fOilVisc;
};

#endif
