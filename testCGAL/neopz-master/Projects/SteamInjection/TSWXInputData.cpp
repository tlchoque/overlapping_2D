//---------------------------------------------------------------------------


#pragma hdrstop

#include "TSWXInputData.h"

//---------------------------------------------------------------------------
double TSWXInputData::getRho_C_Estrela(int i, double So) {
	if (So <= 0. || So >= 1.) {
		return -1.;
	}
	// Determinando o valor da taxa de inje��o da massa de �gua a partir da taxa de inje��o de calor que � assumido constante por partes
	// double massrate = fInjectionData.fTable[i].first/(fWaterInSaturationState.getSaturationStateSpecificEnthalpyToLiquidWater(fSteamTemperature)+fQuality*fWaterInSaturationState.getSaturationStateLatentHeat(fSteamTemperature)-fWaterInSaturationState.getSpecificEnthalpyToLiquidWater(fReservoirData.getReservoirTemperature()));
	double Sw = fWaterInSaturationState.getSaturationWater(i, fInjectionData.getInjectionQuality(), So);
	double val = fReservoirData.Rho_C_Estrela(fWaterInSaturationState, fInjectionData.Temperature(), Sw, So);

	return val;
}


#pragma package(smart_init)
