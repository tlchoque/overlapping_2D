//---------------------------------------------------------------------------

#ifndef TSWXConfinementDataH
#define TSWXConfinementDataH
//---------------------------------------------------------------------------

class TSwxConfinementData
{

	public:

	TSwxConfinementData()
	{
	}

	void SetData(double condtermica, double densidade, double calorespecif);
	TSwxConfinementData(TSwxConfinementData &cpy);

	// M�todo para recuperar o produto das propriedades f�sicas Condutividade_Termica * Densidade * Calor_Especifico
	double getProductOfTheProperties();

	double ThermalConductivity() { return fThermalConductivity; }
	double Density() { return fDensity; }
	double SpecificHeat() { return fSpecificHeat; }

	private :
		/**
	 * fThermalConductivity: Condutividade t�rmica [J/(s m C)]
	 * fDensity: Massa espec�fica [kg/m3]
	 * fSpecificHeat: Calor espec�fico [J/(kg C)]
	 */
	double fThermalConductivity, fDensity, fSpecificHeat;
};

#endif
