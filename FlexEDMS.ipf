#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3	// Use modern global access method and strict wave access.
#pragma IgorVersion = 6.3 // Minimum Igor version required
#pragma version = 1.0-beta.1

// Copyright (c) 2017-2019 Michael C. Heiber
// This source file is part of the FlexEDMS project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The FlexEDMS project can be found on Github at https://github.com/MikeHeiber/FlexEDMS

#include <KBColorizeTraces>

Menu "FlexEDMS"
	"Enter Sample Info", /Q, FEDMS_EnterSampleInfo()
	Submenu "Time of Flight"
		"Load ToF Data Files", /Q, FEDMS_LoadTOFDataFiles()
		"Load ToF Data Folder", /Q, FEDMS_LoadTOFDataFolder(0)
		"Load ToF Intensity Test Folder", /Q, FEDMS_LoadTOFDataFolder(1)
		"Analyze ToF Data", /Q, FEDMS_AnalyzeTOF(0)
		"Analyze ToF Intensity Test Data", /Q, FEDMS_AnalyzeTOF(1)
		"Plot ToF Intensity Test Data", /Q, FEDMS_PlotTOF_IntensityTest()
		"Plot ToF Field Dependence", /Q, FEDMS_PlotTOF_FieldDependences()
		"Plot ToF Temperature Dependence", /Q, FEDMS_PlotTOF_TempDependence()
	End
	Submenu "J-V"
		"Load J-V Data Folder", /Q, FEDMS_LoadJVFolderGUI()
		"Plot J-V Data", /Q, FEDMS_PlotJVGUI()
	End
	Submenu "IPDA"
		"Load J-V Data Folder", /Q, FEDMS_LoadJVFolderGUI()
		"Load Impedance Data Folder", /Q, FEDMS_LoadImpedanceFolderGUI()
		"Calculate IPDA Results", /Q, FEDMS_CalculateIPDAResultsGUI()
		"Calculate and Plot Photocurrent Data", /Q, FEDMS_PlotPhotocurrentDataGUI()
	End
End

Function FEDMS_AnalyzeCapacitanceData(device_name,[show_graphs])
	String device_name
	Variable show_graphs
	String original_data_folder = GetDataFolder(1)
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[Strlen(device_name)-1]
	// Get device info
	SetDataFolder root:FlexEDMS:$(substrate_num)
	NVAR Active_thickness_cm
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance
	Wave/T illuminations = $("illuminations")
	Make/O/D/N=(numpnts(illuminations)) V_sats, V_sat_cors, C_sats, n_sats, n_scs, n_mps, n_ocs, k_mps, k_ocs, mu_mps, V_scs, V_mps, V_ocs
	// Integrate Capacitance vs Bias data for each light intensity
	Variable i
	for(i=0;i<numpnts(illuminations);i+=1)
		FEDMS_CalculateCarrierDensity(device_name,illuminations[i])
		SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance:$("CV_"+illuminations[i])
		Wave carrier_density
		Wave voltage_applied
		Wave voltage_cor
		Duplicate/O voltage_applied $("V_int") $("F_int") $("k_br") $("mu_eff")
		Wave V_int
		Wave F_int
		Wave k_br
		Wave mu_eff
		NVAR V_sat
		NVAR V_sat_cor
		NVAR C_sat
		NVAR n_sat
		NVAR n_sc
		NVAR n_mp
		NVAR n_oc
		V_sats[i] = V_sat
		V_sat_cors[i] = V_sat_cor
		C_sats[i] = C_sat
		n_sats[i] = n_sat
		n_scs[i] = n_sc
		n_mps[i] = n_mp
		n_ocs[i] = n_oc
		// Gather appropriate J-V data
		SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV:$(illuminations[i])
		Wave J_avg
		Wave J_br
		NVAR V_oc
		NVAR V_mp
		V_scs[i] = 0
		V_mps[i] = V_mp
		V_ocs[i] = V_oc
		// Calculate Bias dependent k_br and mu_eff
		Variable j
		for(j=0;j<numpnts(voltage_applied);j+=1)
			V_int[j] = V_oc-voltage_cor[j]
			F_int[j] = V_int[j]/Active_thickness_cm
			k_br[j] = 1e-3*J_br(voltage_applied[j])/(1.602177e-19*Active_thickness_cm*carrier_density[j]^2)
			mu_eff[j] = 1e-3*-J_avg(voltage_applied[j])*Active_thickness_cm/(2*1.60217662e-19*carrier_density[j]*V_int[j])
		endfor
		k_mps[i] = {interp(V_mp,voltage_applied,k_br)}
		k_ocs[i] = {interp(V_oc,voltage_applied,k_br)}
		mu_mps[i] = {interp(V_mp,voltage_applied,mu_eff)}
	endfor
	// Plot calculated charge carrier density graphs
	FEDMS_PlotCarrierDensityVsBias(device_name)
	SetAxis bottom 0.2,*
	SetAxis left 1e15,*
	// Determine if any illuminations should be rejected from further analysis
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance
	Wave rejection_mask
	String graph_name = WinName(0,1)
	ModifyGraph margin(top)=140
	DrawText 0,0,"Select any data sets that should be rejected from further analysis."
	// Add Checkbox for each illumination
	for(i=0;i<numpnts(illuminations);i+=1)
		CheckBox $("checkbox_"+num2str(i)) title=illuminations[i],pos={50,17*i+5},value=rejection_mask[i]
	endfor
	Button button0 title="Continue",proc=FEDMS_Button_RejectIPDA,size={80,25},pos={200,50},fSize=14
	PauseForUser $graph_name
	// Plot final k_br and mu_eff data
	if(show_graphs==1)
		FEDMS_PlotKbrVsBias(device_name)
		FEDMS_PlotMobilityVsBias(device_name)
	endif
	SetDataFolder original_data_folder
End

Function FEDMS_AnalyzeDarkImpedanceCf(device_name,[show_graphs])
	String device_name
	Variable show_graphs
	String original_data_folder = GetDataFolder(1)
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[Strlen(device_name)-1]
	// Gather device information
	SetDataFolder root:FlexEDMS:$(substrate_num)
	NVAR Device_area_cm2
	NVAR Active_thickness_cm
	// Determine the AC series resistance from the dark impedance data
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance:Cf_dark
	Wave frequency
	String wave_list = WaveList("Z_real*",";","")
	wave_list = SortList(wave_list,";",16)
	Variable wave_count = ItemsInList(wave_list)
	Variable x_start = pnt2x(frequency,BinarySearch(frequency,5e5))
	Variable x_end = pnt2x(frequency,BinarySearch(frequency,5e6))
	// Choose which reverse bias to use for calculating the series resistance
	if(wave_count>1)
		Make/D/N=(wave_count)/O R_s_wave
		Make/T/N=(wave_count)/O waves
		Variable i
		for(i=0;i<wave_count;i+=1)
			waves[i] = StringFromList(i,wave_list)
			Wave Z_real = $(StringFromList(i,wave_list))
			R_s_wave[i] = Mean(Z_real,x_start,x_end)
			// Create Graph
			if(i==0)
				Display Z_real vs frequency
			else
				AppendToGraph Z_real vs frequency
			endif
		endfor
		Execute "FEDMS_GraphStyle()"
		SetAxis bottom 100000,*
		SetAxis left 10,400
		ModifyGraph marker=19,msize=1
		ModifyGraph useMrkStrokeRGB=0
		Label left "Z_real (Ohm)"
		Label bottom "Frequency (Hz)"
		Legend/C/N=text0/F=0/A=LB
		TextBox/C/N=text1/F=0/A=MT/X=0.00/Y=6.00 "Dark Impedance"
		// Add listbox and button for user input
		String graph_name = WinName(0,1)
		DoWindow/C $graph_name
		ModifyGraph margin(top)=50
		DrawText 0,0,"\\Z10Select which curve to use to calculate the series resistance."
		ListBox list0 size={90,20*wave_count},listWave=waves,mode=2
		Button button0 title="Go",proc=FEDMS_Button_ChooseListItem,size={60,25},fSize=14
		PauseForUser $graph_name
		NVAR selection_index
		Variable/D/G R_s_ac = R_s_wave[selection_index]
		KillVariables/Z selection_index
		KillWaves/Z R_s_wave waves
	else
		Wave Z_real = $(StringFromList(i,wave_list))
		if(show_graphs==1)
			Display Z_real vs frequency
			Execute "FEDMS_GraphStyle()"
			SetAxis bottom 100000,*
			SetAxis left 10,400
			ModifyGraph marker=19,msize=1
			Label left "Z_real (Ohm)"
			Label bottom "Frequency (Hz)"	
			TextBox/C/N=text1/F=0/A=MT/X=0.00/Y=6.00 "Dark Impedance"
		endif
		Variable/D/G R_s_ac = Mean(Z_real,x_start,x_end)
	endif
	// Calculate the inductance and geometric capacitance from the dark impedance data
	wave_list = WaveList("Z_imag*",";","")
	wave_list = SortList(wave_list,";",16)
	wave_count = ItemsInList(wave_list)
	Make/D/N=(wave_count)/O inductance_wave
	Make/D/N=(wave_count)/O C_g_wave
	Variable start_index, end_index
	for(i=0;i<wave_count;i+=1)
		Wave Z_imag = $(StringFromList(i,wave_list))
		Make/D/N=2/O W_coef
		Make/D/N=2/O W_sigma
		W_coef = {1e-7,1e-9}
		start_index = BinarySearch(frequency,6e5)
		end_index = numpnts(frequency)-1
		Duplicate/O frequency $("fit_Z_imag"+num2str(i+1))
		FuncFit/N/Q/NTHR=0/W=2 FEDMS_FitInductance W_Coef Z_imag[start_index,end_index] /X=frequency /D=$("fit_Z_imag"+num2str(i+1))
		inductance_wave[i] = W_coef[0]
		C_g_wave[i] = W_coef[1]
	endfor
	Variable/G/D L = Mean(inductance_wave)
	Variable/G/D C_g = Mean(C_g_wave)
	// Make a graph to visually check the inductance fits
	if(show_graphs==1)
		for(i=0;i<wave_count;i+=1)
			Wave Z_imag = $(StringFromList(i,wave_list))
			Wave fit = $("fit_Z_imag"+num2str(i+1))
			if(i==0)
				Display Z_imag vs frequency
			else
				AppendToGraph Z_imag vs frequency
			endif
			AppendToGraph fit[start_index,end_index] vs frequency[start_index,end_index]
			ModifyGraph lstyle[2*i+1]=3,rgb[2*i+1]=(0,0,0)
		endfor
		Execute "FEDMS_GraphStyle()"
		ModifyGraph log(left)=0
		SetAxis/A=2 left
		SetAxis bottom 500000,*
		ModifyGraph marker=19,msize=1
		ModifyGraph useMrkStrokeRGB=0
		Label left "Z_imag (Ohm)"
		Label bottom "Frequency (Hz)"	
		Legend/C/N=text0/F=0/A=LT
		TextBox/C/N=text1/F=0/A=MT/X=0.00/Y=6.00 "Dark Impedance"
	endif
	// Calculate the capacitance from the dark impedance data
	wave_list = WaveList("Z_real*",";","")
	wave_list = SortList(wave_list,";",16)
	wave_count = ItemsInList(wave_list)
	Make/T/O/N=(wave_count) waves
	for(i=0;i<wave_count;i+=1)
		Wave Z_real = $(StringFromList(i,wave_list))
		String bias_ending = StringFromList(2,StringFromList(i,wave_list),"_")
		waves[i] = "capacitance_"+bias_ending
		Duplicate/O Z_real $("capacitance_"+bias_ending)
		Wave capacitance = $("capacitance_"+bias_ending)
		Wave Z_imag = $("Z_imag_"+bias_ending)
		capacitance = (-1/(2*PI*frequency))*((Z_imag-2*PI*frequency*L)/((Z_real-R_s_ac)^2+(Z_imag-2*PI*frequency*L)^2))
	endfor
	// Plot the capacitance data
	if(wave_count>1 || show_graphs==1)
		for(i=0;i<wave_count;i+=1)
			bias_ending = StringFromList(2,StringFromList(i,wave_list),"_")
			Wave capacitance = $("capacitance_"+bias_ending)
			if(i==0)
				Display capacitance vs frequency
			else
				AppendToGraph capacitance vs frequency
			endif
		endfor
		Execute "FEDMS_GraphStyle()"
		SetAxis bottom *,*
		ModifyGraph log(left)=0
		SetAxis/A=2 left
		ModifyGraph marker=19,msize=1
		ModifyGraph useMrkStrokeRGB=0
		Label left "Capacitance (F)"
		Label bottom "Frequency (Hz)"	
		Legend/C/N=text0/F=0/A=LB
		TextBox/C/N=text1/F=0/A=MT/X=0.00/Y=6.00 "Dark Capacitance"
	endif
	// Choose which dark capacitance data to use for further analysis
	if(wave_count>1)
		// Add listbox and button for user input
		graph_name = WinName(0,1)
		DoWindow/C $graph_name
		ModifyGraph margin(top)=50
		DrawText 0,0,"\\Z10Select which curve to use to calculate the depleted device capacitance."
		ListBox list0 size={120,20*wave_count},listWave=waves,mode=2
		Button button0 title="Go",proc=FEDMS_Button_ChooseListItem,size={60,25},fSize=14
		PauseForUser $graph_name
		NVAR selection_index
		Variable/D/G C_g = C_g_wave[selection_index]
		Duplicate/O $waves[selection_index] $("capacitance")
		KillVariables/Z selection_index
	else
		Duplicate/O $waves[0] $("capacitance")
		Variable/D/G C_g = C_g_wave[0]
	endif
	// Calculate dielectric constant from the geometric capacitance
	Variable/G/D epsilon = C_g*Active_thickness_cm/(8.854187e-12*1e-2*Device_area_cm2)
	// Clean up
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance:Cf_dark
	KillWaves/Z T_Constraints W_coef W_sigma C_g_wave waves 
	// Restore original working directory
	SetDataFolder original_data_folder
End

Function FEDMS_AnalyzeLightImpedanceCV(device_name,measurement_name,[show_graphs])
	String device_name
	String measurement_name
	Variable show_graphs
	String original_data_folder = GetDataFolder(1)
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[Strlen(device_name)-1]
	// Gather device information
	SetDataFolder root:FlexEDMS:$(substrate_num)
	NVAR Device_area_cm2
	NVAR Active_thickness_cm
	// Gather dark impedance data
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance:Cf_dark
	Wave frequency_dark = $("frequency")
	Wave capacitance_dark = $("capacitance")
	if(!WaveExists(capacitance_dark))
		if(show_graphs==1)
			FEDMS_AnalyzeDarkImpedanceCf(device_name,show_graphs=1)
		else
			FEDMS_AnalyzeDarkImpedanceCf(device_name,show_graphs=0)
		endif
	endif
	Wave frequency_dark = $("frequency")
	Wave capacitance_dark = $("capacitance")
	NVAR L
	NVAR R_s_ac
	// Calculate the dark CV data
	if(DataFolderExists("root:FlexEDMS:"+substrate_num+":"+device_num+":Impedance:CV_dark"))
		SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance:CV_dark
		Wave voltage_applied_dark = $("voltage_applied")
		Duplicate/O voltage_applied_dark $"voltage_cor" $"C_tot"
		Wave Z_real_dark = $("Z_real")
		Wave Z_imag_dark = $("Z_imag")
		NVAR Freq
		Wave C_tot_dark = $("C_tot")
		C_tot_dark = (-1/(2*PI*Freq))*((Z_imag_dark-2*PI*Freq*L)/((Z_real_dark-R_s_ac)^2+(Z_imag_dark-2*PI*Freq*L)^2))
	endif
	// Calculate the capacitance from the light impedance data
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance:$("CV_"+measurement_name)
	Wave voltage_applied
	Duplicate/O voltage_applied C_tot C_mu
	Wave Z_real
	Wave Z_imag
	NVAR Freq
	C_tot = (-1/(2*PI*Freq))*((Z_imag-2*PI*Freq*L)/((Z_real-R_s_ac)^2+(Z_imag-2*PI*Freq*L)^2))
	// Substract the dark capacitance from the light capacitance to determine the chemical capacitance
	FindLevel/P/Q frequency_dark, Freq
	Variable Freq_index = V_LevelX
	Variable C_d = Mean(capacitance_dark,pnt2x(capacitance_dark,Freq_index-2),pnt2x(capacitance_dark,Freq_index+2))
	C_mu = C_tot - C_d
	// Calculate alternate C_mu
	if(DataFolderExists("root:FlexEDMS:"+substrate_num+":"+device_num+":Impedance:CV_dark"))
		Duplicate/O voltage_applied C_mu_alt
		C_mu_alt = C_tot - C_tot_dark
	endif
	// Correct the bias for the series resistance
	//Duplicate/O voltage_applied voltage_cor
	//SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV:$(measurement_name):
	//Wave J_avg
	//Variable i
	//for(i=0;i<numpnts(voltage_applied);i+=1)
		//if(voltage_applied[i]<pnt2x(J_avg,numpnts(J_avg)-1))
			//voltage_cor[i] = voltage_applied[i]-abs(J_avg(voltage_applied[i])*1e-3)*Device_area_cm2*R_s_ac
		//endif
	//endfor
	// Clean up from light analysis
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance:$("CV_"+measurement_name):
	// Restore original working directory
	SetDataFolder original_data_folder
End

Function FEDMS_AnalyzeImpedanceData(device_name,[show_graphs])
	String device_name
	Variable show_graphs
	String original_data_folder = GetDataFolder(1)
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[Strlen(device_name)-1]
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance
	Wave/T illuminations
	Make/N=(numpnts(illuminations))/O intensities
	// Sort illuminations by intensity
	FEDMS_AnalyzeJVIntensity(device_name)
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV:
	Wave/T illuminations_JV = $("illuminations")
	Wave intensities_JV = $("intensities")
	Variable i, j
	for(i=0;i<numpnts(illuminations);i+=1)
		for(j=0;j<numpnts(illuminations_JV);j+=1)
			if(StringMatch(illuminations[i],illuminations_JV[j]))
				intensities[i] = intensities_JV[j]
				break
			endif
		endfor
	endfor
	Sort intensities, illuminations, intensities
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance
	Wave rejection_mask
	if(!WaveExists(rejection_mask))
		Make/D/N=(numpnts(illuminations))/O $("rejection_mask")
		Wave rejection_mask
		rejection_mask = 0
	endif
	// Analyze Impedance data to calculate the chemical capacitance vs corrected bias
	if(show_graphs==1)
		FEDMS_AnalyzeDarkImpedanceCf(device_name,show_graphs=1)
		FEDMS_AnalyzeLightImpedanceCV(device_name,illuminations[0],show_graphs=1)
	else
		FEDMS_AnalyzeDarkImpedanceCf(device_name,show_graphs=0)
		FEDMS_AnalyzeLightImpedanceCV(device_name,illuminations[0],show_graphs=0)
	endif
	for(i=1;i<numpnts(illuminations);i+=1)
		FEDMS_AnalyzeLightImpedanceCV(device_name,illuminations[i])
	endfor
	if(show_graphs==1)
		FEDMS_PlotCapacitanceVsBias(device_name)
	endif
	SetDataFolder original_data_folder
End

Function FEDMS_AnalyzeJV(device_name,measurement_type,measurement_name)
	String device_name
	String measurement_type
	String measurement_name
	String original_data_folder = GetDataFolder(1)
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[Strlen(device_name)-1]
	SetDataFolder root:FlexEDMS:$(substrate_num)
	// Gather device info
	SVAR File_format
	NVAR Device_area_cm2
	NVAR Active_thickness_cm
	// Analyze the designated JV measurement
	SetDataFolder :$(device_num):$(measurement_type):
	NVAR Mismatch_factor
	SetDataFolder $measurement_name
	Wave voltage
	Wave current
	Duplicate/O current J_raw
	if(StringMatch(measurement_name,"*dark*"))
		J_raw = 1000*current/(Device_area_cm2)
	else
		J_raw = 1000*current/(Device_area_cm2*Mismatch_factor)
	endif
	Variable voltage_step = abs(round(10000*(voltage[1]-voltage[0]))/10000)
	Variable voltage_start1
	Variable voltage_start2
	Variable sweep_point1 = 0
	Variable sweep_point2 = 0
	Variable i
	// Separate forward and reverse current sweeps
	if(stringmatch(File_format,"Richter Lab"))
		if(voltage[0]<0)
			for(i=1;i<numpnts(voltage);i+=1)
				if(voltage[i-1]-voltage[i]>0.001)
					sweep_point1 = i-2
					break
				endif
			endfor
			Duplicate/O/R=(0,sweep_point1) J_raw J_forward
			Duplicate/O/R=(sweep_point2,numpnts(voltage)-1) J_raw J_reverse
			voltage_start1 = round(10000*voltage[0])/10000
			voltage_start2 = round(10000*voltage[sweep_point2])/10000
		else // Inverted sweep
			start_point = 0
			for(i=1;i<numpnts(voltage);i+=1)
				if(voltage[i]-voltage[i-1]>0.001)
					start_point = i-2
					break
				endif
			endfor
			Duplicate/O/R=(0,start_point) J_raw J_reverse
			Duplicate/O/R=(start_point+1,numpnts(voltage)-1) J_raw J_forward
			voltage_start2 = round(10000*voltage[0])/10000
			voltage_start1 = round(10000*voltage[sweep_point2])/10000
		endif
	elseif(stringmatch(File_format,"Hersam Lab"))
		// Only works with standard forward+reverse scans that have same voltage range in both directions
		
		sweep_point1 = (numpnts(voltage)/2)-1
		sweep_point2 = (numpnts(voltage)/2)
		Duplicate/O/R=(0,sweep_point1) J_raw J_forward
		Duplicate/O/R=(sweep_point2,numpnts(voltage)-1) J_raw J_reverse
		voltage_start1 = round(10000*voltage[0])/10000
	endif
	SetScale/P x voltage_start1,voltage_step,"", J_forward
	SetScale/P x voltage_start2,(voltage_step*-1),"", J_reverse
	Reverse J_reverse
	Duplicate/O J_forward J_avg
	for(i=0;i<numpnts(J_avg);i+=1)
		J_avg[i] = (J_forward[i] + J_reverse[i])/2
	endfor
	// Determine light J-V characteristics
	Variable index_start, index_end
	if(!StringMatch(measurement_name,"*dark*"))
		FEDMS_CalculatePVJVResults(J_avg)
		WAve W_result
		Variable/G/D J_sc = W_result[0]
		Variable/G/D V_oc = W_result[1]
		Variable/G/D P_max = W_result[2]
		Variable/G/D J_mp = W_result[3]
		Variable/G/D V_mp = W_result[4]
		Variable/G/D FF = W_result[5]
		KillWaves/Z W_result
	// Determine dark J-V characteristics
	else
		Duplicate/O J_avg J_avg_abs
		J_avg_abs = abs(J_avg)
		// Determine shunt resistance
		CurveFit/N/NTHR=0/W=2/Q line J_avg[x2pnt(J_avg,-0.2),x2pnt(J_avg,0.2)] /D
		Wave W_coef
		Variable/G/D R_sh = (1000/(W_coef[1]*Device_area_cm2))
		// Determine series resistance
		WaveStats/Q J_avg
		index_start = numpnts(J_avg)-4
		index_end = numpnts(J_avg)-1
		CurveFit/N/NTHR=0/W=2/Q line J_avg[index_start,index_end] /D
		Variable/G R_s = (1000/(W_coef[1]*Device_area_cm2))
	endif
	// Clean up
	KillWaves/Z J_forward, J_reverse, fit_J_avg, W_coef, W_sigma, W_fitConstants, J_raw
	SetDataFolder $(original_data_folder)
End

Function FEDMS_AnalyzeJVIntensity(device_name)
	String device_name
	String original_data_folder = GetDataFolder(1)
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[Strlen(device_name)-1]
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV
	Wave/T illuminations
	Make/N=(numpnts(illuminations))/O $("I_suns"), J_sat, intensities, J_scs, V_scs, V_mps, V_ocs, PCEs, FFs
	Wave I_suns_wave = $("I_suns")
	FindValue /TEXT="1sun" illuminations
	Variable index_1sun = V_value
	Variable i
	for(i=0;i<numpnts(illuminations);i+=1)
		SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV:$(illuminations[i]):
		FEDMS_AnalyzeJV(device_name,"JV",illuminations[i])
		FEDMS_CalculatePhotocurrent(device_name,illuminations[i])
		NVAR V_0
		Wave photocurrent
		J_sat[i] = photocurrent(V_0-2)
	endfor
	I_suns_wave = J_sat/J_sat[index_1sun]
	intensities = 100*I_suns_wave
	Sort I_suns_wave, illuminations, intensities, I_suns_wave, J_sat
	for(i=0;i<numpnts(illuminations);i+=1)
		SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV:$(illuminations[i]):
		Variable/G I_suns = I_suns_wave[i]
		NVAR J_sc
		NVAR V_mp
		NVAR V_oc
		NVAR P_max
		NVAR FF
		J_scs[i] = J_sc
		V_scs[i] = 0
		V_mps[i] = V_mp
		V_ocs[i] = V_oc
		PCEs[i] = 100*P_max/intensities[i]
		FFs[i] = FF
	endfor
	// Calculate slope of ln(I) vs V_oc
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV
	Duplicate/O intensities ln_I
	ln_I = ln(intensities)
	CurveFit/N/Q/NTHR=0/W=2 line V_ocs[0,index_1sun] /X=ln_I /D
	Wave W_coef
	Variable/D/G ideality_factor = W_coef[1]/(8.6173e-5*298)
	// Clean up
	KillWaves/Z J_sat ln_I W_coef $"W_sigma" $"fitX_V_ocs" $"fit_V_ocs"
	SetDataFolder original_data_folder
End

Function FEDMS_AnalyzePhotocurrent(device_name,measurement_name)
	String device_name
	String measurement_name
	String original_data_folder = GetDataFolder(1)
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[strlen(device_name)-1]
	// Check that the device folder and measurement folder exist
	if(!DataFolderExists("root:FlexEDMS:"+substrate_num+":"+device_num+":JV:'"+measurement_name+"'"))
		Print "Error! Device name or measurement not found. Check that the sample and data have been correctly loaded."
		return NaN
	endif
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV:$(measurement_name)
	Variable/G V_status = 1
	// Fit photocurrent in saturation regime
	Wave photocurrent
	Wave voltage_effective
	Display photocurrent vs voltage_effective
	Execute "FEDMS_GraphStyle()"
	ModifyGraph expand=1.2, log(bottom)=1, rgb[0]=(32768,40777,65535), useMrkStrokeRGB[0]=0
	ModifyGraph log(left)=0
	SetAxis left 0,*
	SetAxis bottom 0.04,*
	TextBox/C/N=text0/F=0/Z=1/A=LT device_name+", "+measurement_name
	String trace_name = StringFromList(0,TraceNameList("",";",1))
	// Plot the photocurrent fit from the previous analysis if it exists
	Wave fit_photocurrent
	if(WaveExists(fit_photocurrent))
		AppendToGraph fit_photocurrent vs voltage_effective
		ModifyGraph lstyle(fit_photocurrent)=3,rgb(fit_photocurrent)=(0,0,0)
		NVAR A_position
		NVAR B_position
		Cursor/P A $trace_name A_position
		Cursor/P B $trace_name B_position
	else
		Cursor A $trace_name 0.15
		Cursor/P B $trace_name numpnts(photocurrent)-1
	endif
	ShowInfo
	String graph_name = WinName(0,1)
	DoWindow/C $graph_name
	ModifyGraph margin(left)=35
	ModifyGraph margin(top)=25
	DrawText 0,0,"\\Z10Move cursors to select fit range and click Fit."
	Button button0 title="Fit",proc=FEDMS_Button_FitPhotocurrent,size={60,25},fSize=14
	Button button1 title="Done",proc=FEDMS_Button_Done,size={60,25},fSize=14
	Button button2 title="Cancel",proc=FEDMS_Button_Cancel,size={60,25},fSize=14
	PauseForUser $graph_name
	// Check if user has cancelled the fit
	if(V_status==-1)
		KillVariables/Z V_status
		SetDataFolder original_data_folder
		Return NaN
	endif
	// Get final photocurrent fit
	Wave fit_photocurrent
	// Determine V_sat_onset characteristics
	Variable/G V_sat_onset
	Variable/G V_eff_onset
	Variable/G J_sat_onset
	NVAR V_0
	Variable i
	for(i=0;i<numpnts(photocurrent);i+=1)
		if(abs(mean(photocurrent,pnt2x(photocurrent,i-3),pnt2x(photocurrent,i+3))-fit_photocurrent[i])<0.005*abs(fit_photocurrent[i]))
			V_sat_onset = pnt2x(photocurrent,i)
			V_eff_onset = V_0-V_sat_onset
			J_sat_onset = photocurrent[i]
			break
		endif
	endfor
	// Calculate Recombination Current
	Wave J_avg
	Duplicate/O photocurrent, J_gen
	Duplicate/O J_avg, J_br
	J_gen = fit_photocurrent
	Reverse J_gen
	// Mirror J_gen for voltages above V_0
	for(i=0;i<numpnts(J_gen);i+=1)
		if(pnt2x(J_gen,i)>V_0)
			J_gen[i] = J_gen(V_0-(pnt2x(J_gen,i)-V_0))
			if(numtype(J_gen[i])!=0)
				J_gen[i] = J_gen[i-1]
			endif
		endif
	endfor
	J_br = J_gen + J_avg
	// Clean up
	KillVariables/Z V_status $"V_sat"
	KillWaves/Z T_Constraints $("W_sigma") $"J_ideal" $"W_fitConstants"
	SetDataFolder original_data_folder
	return 1
End

Function FEDMS_AnalyzePhotocurrentData(device_name)
	String device_name
	String measurement_name
	String original_data_folder = GetDataFolder(1)
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[strlen(device_name)-1]
	// Check that the device folder and measurement folder exist
	if(!DataFolderExists("root:FlexEDMS:"+substrate_num+":"+device_num))
		Print "Error! Device name not found. Check that the sample and data have been correctly loaded."
		return NaN
	endif
	// Show all photocurrent data
	// Show fits if analyzed previously
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV:$("1sun")
	Wave/Z fit_photocurrent
	if(WaveExists(fit_photocurrent))
		FEDMS_PlotPhotocurrentData(device_name,show_fits=1)
	else
		FEDMS_PlotPhotocurrentData(device_name,show_fits=0)
	endif
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV
	Wave/T illuminations
	Wave intensities
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV
	KillWaves/Z $"J_br_intensity"  $"J_gen_intensity"
	Make/D/O/N=(numpnts(illuminations)) V_sat_onsets, V_eff_onsets, J_sat_onsets
	Variable i
	for(i=0;i<numpnts(illuminations);i+=1)
		Variable status = FEDMS_AnalyzePhotocurrent(device_name,illuminations[i])
		if(numtype(status)==2)
			return NaN
		endif
		SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV:$(illuminations[i])
		NVAR V_sat_onset
		NVAR V_eff_onset
		NVAR J_sat_onset
		V_sat_onsets[i] = V_sat_onset
		V_eff_onsets[i] = V_eff_onset
		J_sat_onsets[i] = J_sat_onset
	endfor
	SetDataFolder original_data_folder
End

Function FEDMS_AnalyzeTOF(isIntensityTest)
	Variable isIntensityTest
	String original_folder = GetDataFolder(1)
	String sample_name = FEDMS_ChooseTOFSample()
	if(StringMatch(sample_name,""))
		return NaN
	endif
	String carrier_type = FEDMS_ChooseTOFCarrierType(sample_name)
	if(StringMatch(carrier_type,""))
		return NaN
	endif
	SetDataFolder root:FlexEDMS:$(sample_name):$("Time of Flight"):$(carrier_type)
	if(isIntensityTest)
		SetDataFolder root:FlexEDMS:$(sample_name):$("Time of Flight"):$(carrier_type):$("intensity_test")
	endif
	// Build the temperature list
	String temperature_option
	DFREF dfr1 = GetDataFolderDFR()
	Variable N_temps = CountObjectsDFR(dfr1,4)
	String temperature_list = "all"
	String temperature_name
	Variable i
	for(i=0;i<N_temps;i+=1)
		temperature_name = GetIndexedObjNameDFR(dfr1,4,i)
		if(StringMatch(temperature_name,"!intensity_test"))
			temperature_list = AddListItem(temperature_name,temperature_list)
		endif
	endfor
	temperature_list = SortList(temperature_list,";",16)
	// Prompt user to choose temperature
	Prompt temperature_option, "Choose a temperature option:" popup, temperature_list
	DoPrompt "Make Selection", temperature_option
	if(V_flag==1)
		SetDataFolder original_folder
		return NaN
	endif
	// Create needed data waves and variables
	NVAR Thickness_m = root:FlexEDMS:$(sample_name):Thickness_m
	for(i=0;i<N_temps;i+=1)
		Variable measurement_count = 0
		temperature_name = StringFromList(i+1,temperature_list)
		if(StringMatch(temperature_option,"all") || StringMatch(temperature_name,temperature_option))
			SetDataFolder :$(temperature_name)
			Wave/Z voltage_V
			if(!WaveExists(voltage_V))
				Make/N=1/D $"voltage_V"
				Wave voltage_V
			endif
			Wave/Z electric_field
			if(!WaveExists(electric_field))
				Make/N=1/D $"electric_field"
				Wave electric_field
			endif
			Wave/Z transit_time_geo_s
			if(!WaveExists(transit_time_geo_s))
				Make/N=1/D $"transit_time_geo_s"
				Wave transit_time_geo_s
			endif
			Wave/Z transit_time_half_s
			if(!WaveExists(transit_time_half_s))
				Make/N=1/D $"transit_time_half_s"
				Wave transit_time_half_s
			endif
			Wave/Z dispersion
			if(!WaveExists(dispersion))
				Make/N=1/D $"dispersion"
				Wave dispersion
			endif
			Wave/Z mobility_geo
			if(!WaveExists(mobility_geo))
				Make/N=1/D $"mobility_geo"
				Wave mobility_geo
			endif
			Wave/Z rejection_mask
			if(!WaveExists(rejection_mask))
				Make/N=1/D $"rejection_mask"
				Wave rejection_mask
			endif
			Wave/Z intensity
			if(isIntensityTest && !WaveExists(intensity))
				Make/N=1/D $"intensity"
				Wave intensity
			endif
			// Determine number of biases and build bias list
			DFREF dfr2 = GetDataFolderDFR()
			Variable N_biases = CountObjectsDFR(dfr2,4)
			String bias_list = ""
			Variable j
			for(j=0;j<N_biases;j+=1)
				bias_list = AddListItem(GetIndexedObjNameDFR(dfr2,4,j),bias_list)
			endfor
			if(StringMatch(StringFromList(0,bias_list),"-*"))
				bias_list = SortList(bias_list,";",17)
			else
				bias_list = SortList(bias_list,";",16)
			endif
			// Loop through all biases for the given temperature
			for(j=0;j<N_biases;j+=1)
				String bias_name = StringFromList(j,bias_list)
				SetDataFolder :$(bias_name)
				// If analyzing intensity test data, build the intensity list
				if(isIntensityTest)
					DFREF dfr3 = GetDataFolderDFR()
					Variable N_intensities = CountObjectsDFR(dfr3,4)
					String intensity_list = ""
					Variable k
					for(k=0;k<N_intensities;k+=1)
						intensity_list = AddListItem(GetIndexedObjNameDFR(dfr3,4,k),intensity_list)
					endfor
					intensity_list = SortList(intensity_list,";",16)
					for(k=0;k<N_intensities;k+=1)
						String intensity_name = StringFromList(k,intensity_list)
						SetDataFolder :$(intensity_name)
						Variable/G V_status = 1
						Variable/G Measurement_index = measurement_count
						intensity[Measurement_index] = {str2num(RemoveEnding(intensity_name,"%"))}
						voltage_V[Measurement_index] = {str2num(RemoveEnding(bias_name,"V"))}
						electric_field[Measurement_index] = {abs(voltage_V[Measurement_index])/Thickness_m}
						FEDMS_AnalyzeTOFData(GetDataFolder(1),sample_name,temperature_name,carrier_type,bias_name)
						if(V_status==-1)
							KillVariables V_status
							SetDataFolder original_folder
							return NaN
						endif
						FEDMS_AnalyzeTOF_TempDependence(sample_name,carrier_type)
						measurement_count+=1
						KillVariables V_status
						SetDataFolder ::
					endfor
				else
					Variable/G V_status = 1
					Variable/G Measurement_index = measurement_count
					voltage_V[Measurement_index] = {str2num(RemoveEnding(bias_name,"V"))}
					electric_field[Measurement_index] = {abs(voltage_V[Measurement_index])/Thickness_m}
					FEDMS_AnalyzeTOFData(GetDataFolder(1),sample_name,temperature_name,carrier_type,bias_name)
					if(V_status==-1)
						KillVariables V_status
						SetDataFolder original_folder
						return NaN
					endif
					FEDMS_AnalyzeTOF_TempDependence(sample_name,carrier_type)
					measurement_count+=1
					KillVariables V_status
				endif
				SetDataFolder ::
			endfor
			SetDataFolder ::
		endif
	endfor
	SetDataFolder original_folder
End

Function FEDMS_AnalyzeTOFData(data_path,sample_name,temperature_name,carrier_type,bias_name)
	String data_path
	String sample_name
	String temperature_name
	String carrier_type
	String bias_name
	String original_folder = GetDataFolder(1)
	Variable isIntensityTest = StringMatch(StringFromList(5,data_path,":"),"intensity_test")
	NVAR Measurement_index
	if(isIntensityTest)
		SetDataFolder root:FlexEDMS:$(sample_name):$("Time of Flight"):$(carrier_type):$("intensity_test"):$(temperature_name)
	else
		SetDataFolder root:FlexEDMS:$(sample_name):$("Time of Flight"):$(carrier_type):$(temperature_name)
	endif
	Wave rejection_mask
	SetDataFolder data_path
	FEDMS_PlotTOFData(data_path,sample_name,carrier_type,temperature_name,bias_name)
	ModifyGraph expand=1.3
	String trace_name = StringFromList(0,TraceNameList("",";",1))
	Wave times = XWaveRefFromTrace("",trace_name)
	Wave fit_pre_transit
	Wave fit_post_transit
	if(WaveExists(fit_pre_transit))
		AppendToGraph fit_pre_transit vs times
		ModifyGraph lstyle(fit_pre_transit)=3,rgb(fit_pre_transit)=(0,0,0)
		NVAR A_position
		Cursor A $trace_name A_position
	else
		Cursor A $trace_name BinarySearch(times,2e-6)
	endif
	if(WaveExists(fit_post_transit))
		AppendToGraph fit_post_transit vs times
		ModifyGraph lstyle(fit_post_transit)=3,rgb(fit_post_transit)=(0,0,0)
		NVAR B_position
		Cursor B $trace_name B_position
	else
		Cursor B $trace_name BinarySearch(times,2e-5)
	endif
	ShowInfo
	String graph_name = WinName(0,1)
	DoWindow/C $graph_name
	ModifyGraph margin(left)=35
	ModifyGraph margin(top)=25
	if(rejection_mask[Measurement_index])
		TextBox/C/N=text1/F=0/Z=1/A=LB "Previously rejected"
	endif
	DrawText 0,0,"\\Z10Move cursors to pre- and post-transit regions and click Analyze."
	Button button0 title="Analyze",proc=FEDMS_ButtonProc_FitTOF,size={60,25},fSize=14
	Button button1 title="Done",proc=FEDMS_ButtonProc_DoneTOF,size={60,25},fSize=14
	Button button2 title="Cancel",proc=FEDMS_ButtonProc_CancelTOF,size={60,25},fSize=14
	Button button3 title="Reject",proc=FEDMS_ButtonProc_RejectToF,size={60,25},fSize=14
	PauseForUser $graph_name
	SetDataFolder original_folder
End

Function FEDMS_AnalyzeTOF_TempDependence(sample_name,carrier_type)
	String sample_name
	String carrier_type
	String original_folder = GetDataFolder(1)
	SetDataFolder root:FlexEDMS:$(sample_name):$("Time of Flight"):$(carrier_type)
	// Build the temperature list
	String temperature_option
	DFREF dfr1 = GetDataFolderDFR()
	Variable N_temps = CountObjectsDFR(dfr1,4)
	Make/N=(N_temps)/D/O temperature_K, mobility_avg, mobility_stdev, dispersion_avg, dispersion_stdev
	String temperature_list = ""
	Variable i
	for(i=0;i<N_temps;i+=1)
		String temperature_name = GetIndexedObjNameDFR(dfr1,4,i)
		if(StringMatch(temperature_name,"!intensity_test"))
			temperature_list = AddListItem(temperature_name,temperature_list)
		endif
	endfor
	temperature_list = SortList(temperature_list,";",16)
	Make/D/N=5/O mobility_temp
	for(i=0;i<N_temps;i+=1)
		temperature_name = StringFromList(i,temperature_list)
		SetDataFolder :$(temperature_name)
		temperature_K[i] = str2num(RemoveEnding(temperature_name,"K"))
		Wave mobility_geo
		Wave dispersion
		if(!WaveExists(mobility_geo))
			mobility_avg[i] = NaN
			mobility_stdev[i] = NaN
			dispersion_avg[i] = NaN
			dispersion_stdev[i] = NaN
			continue
		endif
		//Duplicate/O dispersion dispersion_temp
		//Variable j
		//for(j=0;j<numpnts(mobility_temp);j+=1)
		//	Variable k
		//	Variable current_val = 1
		//	Variable index = 0
		//	for(k=0;k<numpnts(dispersion_temp);k+=1)
		//		if(dispersion_temp[k]>0 && dispersion_temp[k]<current_val)
		//			current_val = dispersion_temp[k]
		//			index = k
		//		endif
		//	endfor
		//	mobility_temp[j] = mobility_geo[index]
		//	dispersion_temp[index] = 1
		//endfor
		//WaveStats/Q mobility_temp
		WaveStats/Q mobility_geo
		mobility_avg[i] = V_avg
		mobility_stdev[i] = V_sdev
		WaveStats/Q dispersion
		dispersion_avg[i] = V_avg
		dispersion_stdev[i] = V_sdev
		//KillWaves dispersion_temp
		SetDataFolder ::
	endfor
	KillWaves mobility_temp
	SetDataFolder original_folder
End

Function FEDMS_Button_Cancel(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			NVAR V_status
			V_status = -1
			KillWindow $WinName(0,1)
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function FEDMS_Button_ChooseListItem(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			ControlInfo list0
			Variable/D/G selection_index = V_value
			if(V_value>=0)
				KillWindow $WinName(0,1)
			endif
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function FEDMS_Button_Done(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			KillWindow $WinName(0,1)
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function FEDMS_Button_FitPhotocurrent(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			String trace_list = TraceNameList("",";",1)
			String trace_name = StringFromList(0,trace_list)
			Wave trace_wave = $trace_name
			Wave voltage_effective
			Duplicate/O trace_wave fit_photocurrent
			Variable/G A_position = pcsr(A)
			Variable/G B_position = pcsr(B)
			Make/O/T/N=3 T_Constraints
			T_Constraints[0] = {"K0 > 0","K1 > 0","K2 > 0"}
			K0=0.5*trace_wave[numpnts(trace_wave)-1];K1=1;K2=0.4
			CurveFit/N/NTHR=0/W=2/Q/G Power, trace_wave[pcsr(A),pcsr(B)] /X=voltage_effective /D=fit_photocurrent /C=T_Constraints
			Wave W_coef
			fit_photocurrent = W_coef[0]+W_coef[1]*voltage_effective^W_coef[2]
			if(ItemsInList(trace_list)==1)
				AppendToGraph fit_photocurrent vs voltage_effective
				ModifyGraph lstyle[1]=3,rgb[1]=(0,0,0)
			endif
			TextBox/C/N=text1/F=0/Z=1/A=RB "J0 = "+num2str(W_coef[0])+"\rJ1 = "+num2str(W_coef[1])+"\rp = "+num2str(W_coef[2])
			// Clean up from fitting
			KillWaves/Z W_coef $"W_sigma" T_Constraints
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function FEDMS_Button_FitTOF(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			String trace_list = TraceNameList("",";",1)
			String trace_name = StringFromList(0,trace_list)
			Wave trace_wave = $trace_name
			Wave times = XWaveRefFromTrace("",trace_name)
			Duplicate/O trace_wave fit_pre_transit fit_post_transit fit_pre_transit_half
			Duplicate/O trace_wave log_trace_wave log_times
			log_trace_wave = log(trace_wave)
			log_times = log(times)
			Variable/G A_position = pcsr(A)
			Variable/G B_position = pcsr(B)
			CurveFit/N/Q line log_trace_wave[pcsr(A)-5,pcsr(A)+5] /X=log_times /D
			Wave W_coef
			fit_pre_transit = 10^(W_coef[0]+W_coef[1]*log_times)
			Duplicate/O W_coef W_coef_pre
			if(ItemsInList(trace_list)==1)
				AppendToGraph fit_pre_transit vs times
				ModifyGraph lstyle(fit_pre_transit)=3,rgb(fit_pre_transit)=(0,0,0)
			endif
			CurveFit/N/Q line log_trace_wave[pcsr(B)-5,pcsr(B)+5] /X=log_times /D
			fit_post_transit = 10^(W_coef[0]+W_coef[1]*log_times)
			Duplicate/O W_coef W_coef_post
			if(ItemsInList(trace_list)==1)
				AppendToGraph fit_post_transit vs times
				ModifyGraph lstyle(fit_post_transit)=3,rgb(fit_post_transit)=(0,0,0)
			endif
			fit_pre_transit_half = 0.5*fit_pre_transit
			if(ItemsInList(trace_list)<4)
				AppendToGraph fit_pre_transit_half vs times
				ModifyGraph lstyle(fit_pre_transit_half)=3,rgb(fit_pre_transit_half)=(0,0,0)
			endif
			// Extract transit times
			Variable transit_time_geo
			Variable transit_time_half
			Variable dispersion_val
			Duplicate/O fit_pre_transit difference
			difference = fit_pre_transit-fit_post_transit
			FindLevel/Q difference, 0
			transit_time_geo = times(V_LevelX)
			difference = fit_pre_transit_half-fit_post_transit
			FindLevel/Q difference, 0
			transit_time_half = times(V_LevelX)
			dispersion_val = (transit_time_half-transit_time_geo)/transit_time_half
			TextBox/C/N=text1/F=0/Z=1/A=LB "t\\Bgeo\\M = "+num2str(transit_time_geo)+"\rt\\B1/2\\M = "+num2str(transit_time_half)+"\rw = "+num2str(dispersion_val)
			// Clean up from fitting
			KillWaves difference $("W_coef") $("W_sigma") log_trace_wave log_times
			// Enter results into table
			String measurement_folder = GetDataFolder(1)
			Variable isIntensityTest = StringMatch(StringFromList(5,measurement_folder,":"),"intensity_test")
			NVAR Measurement_index
			if(isIntensityTest)
				SetDataFolder :::::::
			else
				SetDataFolder :::::
			endif
			NVAR Thickness_m
			SetDataFolder measurement_folder
			if(isIntensityTest)
				SetDataFolder :::
			else
				SetDataFolder ::
			endif
			Wave voltage_V
			Wave transit_time_geo_s
			Wave transit_time_half_s
			Wave dispersion
			Wave mobility_geo
			Wave rejection_mask
			transit_time_geo_s[Measurement_index] = {transit_time_geo}
			transit_time_half_s[Measurement_index] = {transit_time_half}
			dispersion[Measurement_index] = {dispersion_val}
			mobility_geo[Measurement_index] = {(Thickness_m*1e2)^2/(abs(voltage_V[Measurement_index])*transit_time_geo)}
			rejection_mask[Measurement_index] = {0}
			SetDataFolder measurement_folder
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function FEDMS_Button_RejectIPDA(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			Wave rejection_mask
			// Get status of each checkbox
			Variable i
			for(i=0;i<numpnts(rejection_mask);i+=1)
				ControlInfo $("checkbox_"+num2str(i))
				rejection_mask[i] = V_value
			endfor
			KillWindow $WinName(0,1)
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function FEDMS_Button_RejectTOF(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			String measurement_folder = GetDataFolder(1)
			Variable isIntensityTest = StringMatch(StringFromList(5,measurement_folder,":"),"intensity_test")
			NVAR Measurement_index
			if(isIntensityTest)
				SetDataFolder :::
			else
				SetDataFolder ::
			endif
			Wave rejection_mask
			Wave mobility_geo
			Wave dispersion
			rejection_mask[Measurement_index] = {1}
			mobility_geo[Measurement_index] = {NaN}
			dispersion[Measurement_index] = {NaN}
			KillWindow $WinName(0,1)
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function FEDMS_CalculateCarrierDensity(device_name,measurement_name)
	String device_name
	String measurement_name
	String original_data_folder = GetDataFolder(1)
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[Strlen(device_name)-1]
	// Gather device information
	SetDataFolder root:FlexEDMS:$(substrate_num):
	NVAR Device_area_cm2
	NVAR Active_thickness_cm
	// Get J-V curve data
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV:$(measurement_name):
	NVAR V_oc
	NVAR V_mp
	NVAR V_sat_JV = $("V_sat_onset")
	Wave J_avg
	Wave fit_photocurrent
	Wave voltage_effective
	// Get impedance data
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance:$("CV_"+measurement_name):
	Wave C_mu
	//Wave voltage_cor
	Wave voltage_applied
	Duplicate/O voltage_applied carrier_density
	// Determine final V_sat for the impedance data
	Variable/G V_sat
	//Variable/G V_sat_cor
	Variable/G C_sat
	Variable sat_index
	Variable i
	for(i=0;i<numpnts(voltage_applied);i+=1)
		if(voltage_applied[i]>V_sat_JV)
			V_sat = voltage_applied[i]
			//V_sat_cor = voltage_cor[i]
			C_sat = C_mu[i]
			sat_index = i
			break
		endif
	endfor
	Variable/G dmu_sat = 0
	Variable/G mu_sat = 1e-4
	Variable/G dJ_sat = (fit_photocurrent(V_sat+0.05)-fit_photocurrent(V_sat-0.05))/0.1
	Variable/G J_sat = fit_photocurrent(V_sat)
	Variable/G n_sat
	Variable N_pass = 0
	do
		if(N_pass==0)
			Make/N=0/D/O n_sat_wave, mu_sat_wave, dmu_sat_wave
		endif
		//n_sat = (C_sat/(1.602177e-19*Device_area_cm2*Active_thickness_cm))*(1/(V_oc-V_sat_cor)+dJ_sat/J_sat-dmu_sat/mu_sat)^(-1)
		n_sat = (C_sat/(1.602177e-19*Device_area_cm2*Active_thickness_cm))*(1/(V_oc-V_sat)+dJ_sat/J_sat-dmu_sat/mu_sat)^(-1)
		if(n_sat<0)
			n_sat = 0
		endif
		//for(i=0;i<numpnts(voltage_cor);i+=1)
		//	if(voltage_cor[i]<V_sat_cor)
		//		carrier_density[i] = n_sat - AreaXY(voltage_cor,C_mu,voltage_cor[i],V_sat_cor)/(1.602177e-19*Device_area_cm2*Active_thickness_cm)
		//	elseif(voltage_cor[i]>=V_sat_cor)
		//		carrier_density[i] = n_sat + AreaXY(voltage_cor,C_mu,V_sat_cor,voltage_cor[i])/(1.602177e-19*Device_area_cm2*Active_thickness_cm)
		//	endif
		//endfor
		for(i=0;i<numpnts(voltage_applied);i+=1)
			if(voltage_applied[i]<V_sat)
				carrier_density[i] = n_sat - AreaXY(voltage_applied,C_mu,voltage_applied[i],V_sat)/(1.602177e-19*Device_area_cm2*Active_thickness_cm)
			elseif(voltage_applied[i]>=V_sat)
				carrier_density[i] = n_sat + AreaXY(voltage_applied,C_mu,V_sat,voltage_applied[i])/(1.602177e-19*Device_area_cm2*Active_thickness_cm)
			endif
		endfor
		dmu_sat_wave[N_pass] = {dmu_sat}
		n_sat_wave[N_pass] = {n_sat}
		mu_sat = 1e-3*-mean(J_avg,V_sat-0.05,V_sat+0.05)*Active_thickness_cm/(2*1.60217662e-19*n_sat*(V_oc-V_sat))
		mu_sat_wave[N_pass] = {mu_sat}
		// Check for convergence
		if((N_pass>3) && abs(n_sat_wave[N_pass]-n_sat_wave[N_pass-1])>abs(n_sat_wave[N_pass-1]-n_sat_wave[N_pass-2]))
			break
		endif
		if((N_pass>3) && abs(n_sat_wave[N_pass]-n_sat_wave[N_pass-1])<abs(0.0005*n_sat_wave[N_pass]))
			break
		endif
		if(N_pass>500)
			Print "Error! Maximum iterations reached for calculating the carrier density of device "+device_name+" measured at "+measurement_name+"."
			break
		endif
		// Calculate updated mobility derivative for next pass
		Variable mobility_bias1, mobility_bias2
		if(sat_index==0)
			mobility_bias1 = 1e-3*-mean(J_avg,voltage_applied[sat_index]-0.05,voltage_applied[sat_index]+0.05)*Active_thickness_cm/(2*1.60217662e-19*carrier_density[sat_index]*(V_oc-voltage_applied[sat_index]))
			mobility_bias2 = 1e-3*-mean(J_avg,voltage_applied[sat_index+1]-0.05,voltage_applied[sat_index+1]+0.05)*Active_thickness_cm/(2*1.60217662e-19*carrier_density[sat_index+1]*(V_oc-voltage_applied[sat_index+1]))
			dmu_sat = ((mobility_bias2-mobility_bias1)/(voltage_applied[sat_index+1]-voltage_applied[sat_index])+dmu_sat_wave[N_pass])/2
		else
			mobility_bias1 = 1e-3*-mean(J_avg,voltage_applied[sat_index-1]-0.05,voltage_applied[sat_index-1]+0.05)*Active_thickness_cm/(2*1.60217662e-19*carrier_density[sat_index-1]*(V_oc-voltage_applied[sat_index-1]))
			mobility_bias2 = 1e-3*-mean(J_avg,voltage_applied[sat_index+1]-0.05,voltage_applied[sat_index+1]+0.05)*Active_thickness_cm/(2*1.60217662e-19*carrier_density[sat_index+1]*(V_oc-voltage_applied[sat_index+1]))
			dmu_sat = ((mobility_bias2-mobility_bias1)/(voltage_applied[sat_index+1]-voltage_applied[sat_index-1])+dmu_sat_wave[N_pass])/2
		endif
		N_pass += 1
	while(1)
	// Check for a negative n_sat
	if(n_sat<0)
		Print "Error! "+device_name+" measured at "+measurement_name+" resulted in a negative n_sat value."
	endif
	// Save final characteristics
	Variable/G n_sc = interp(0,voltage_applied,carrier_density)
	Variable/G n_mp = interp(V_mp,voltage_applied,carrier_density)
	Variable/G n_oc = interp(V_oc,voltage_applied,carrier_density)
	// Cleanup
	KillWaves/Z $"n_oc_wave"
	// Restore original working directory
	SetDataFolder original_data_folder
End

Function FEDMS_CalculateIPDAFieldDepends(device_name)
	String device_name
	String original_data_folder = GetDataFolder(1)
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[Strlen(device_name)-1]
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance
	Wave/T illuminations
	Wave rejection_mask
	// Calculate electric field dependencies
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):IPDA_Results
	Wave I_suns
	Make/N=21/D/O n_field
	Make/N=(0,numpnts(n_field))/D/O field_sqrt, kbr_vs_field, mu_vs_field
	Variable i
	for(i=0;i<numpnts(n_field);i+=1)
		Variable j
		Variable index = 0
		for(j=0;j<numpnts(illuminations);j+=1)
			if(rejection_mask[j]==1)
				continue
			endif
			SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance:$("CV_"+illuminations[j])
			Wave voltage_applied
			Wave carrier_density
			Wave F_int
			WAve k_br
			Wave mu_eff
			n_field[i] = 10^(0.1+0.1*i)*1e15
			field_sqrt[index][i] = {sqrt(interp(n_field[i],carrier_density,F_int))}
			Variable field_sqrt_limit = sqrt(interp(0,voltage_applied,F_int))
			if(field_sqrt[index][i] < field_sqrt_limit)
				kbr_vs_field[index][i] = {interp(n_field[i],carrier_density,k_br)}
				mu_vs_field[index][i] = {interp(n_field[i],carrier_density,mu_eff)}
			else
				kbr_vs_field[index][i] = {NaN}
				mu_vs_field[index][i] = {NaN}
			endif
			index += 1
		endfor
	endfor
	// Fit field dependent mobility data
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):IPDA_Results
	Make/D/O/N=0 mu_0s, betas
	for(i=0;i<numpnts(n_field);i+=1)
		Make/D/O/N=(dimsize(mu_vs_field,0)) field_data, mobility_data
		Variable index_start = -1
		Variable index_end = -1
		for(j=0;j<dimsize(mu_vs_field,0);j+=1)
			field_data[j] = field_sqrt[j][i]
			if(index_start<0 && numtype(mu_vs_field[j][i])==0 && mu_vs_field[j][i]>0)
				index_start = j
			endif
			if(index_start>=0 && numtype(mu_vs_field[j][i])==0 && mu_vs_field[j][i]>0)
				index_end = j
			endif
			if(numtype(mu_vs_field[j][i])!=0 || mu_vs_field[j][i]<0)
				mobility_data[j] = 0
			else
				mobility_data[j] = mu_vs_field[j][i]
			endif
		endfor
		if(index_end-index_start>1)
			Make/D/N=2/O W_coef
			W_coef[0] = {mobility_data[index_start],0}
			FuncFit/N/NTHR=0/W=2/Q/G FEDMS_FitMobilityField W_coef mobility_data[index_start,index_end] /X=field_data /D
			mu_0s[i] = {W_coef[0]}
			betas[i] = {W_coef[1]}
		else
			mu_0s[i] = {NaN}
			betas[i] = {NaN}
		endif
	endfor
	// Clean up
	KillWaves/Z $("fit_mobility_data") mobility_data field_data
	SetDataFolder original_data_folder
End

Function FEDMS_CalculateIPDAResults(device_name,[hide_table])
	String device_name
	Variable hide_table
	String original_data_folder = GetDataFolder(1)	
	// Perform requisite data analysis
	FEDMS_AnalyzeJVIntensity(device_name)
	FEDMS_AnalyzePhotocurrentData(device_name)
	//FEDMS_PlotPhotocurrentData(device_name,show_fits=1)
	FEDMS_AnalyzeImpedanceData(device_name,show_graphs=0)
	//FEDMS_PlotCapacitanceVsBias(device_name)
	FEDMS_AnalyzeCapacitanceData(device_name,show_graphs=0)
	// Continue with final results calculation
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[Strlen(device_name)-1]
	SetDataFolder root:FlexEDMS:$(substrate_num)
	// Gather device information
	NVAR Device_area_cm2
	NVAR Active_thickness_cm
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance:Cf_dark
	NVAR epsilon
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance
	Wave/T illuminations
	Wave rejection_mask
	Wave intensities
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num)
	NewDataFolder/O/S IPDA_Results
	Make/D/O/N=0 I_suns, n_mps, theta, theta_alpha, gamma_dir, FFs, J_scs, V_ocs, PCEs, tau_ex, tau_rec
	Make/D/O/N=0 G_sc, G_mp, G_oc, V_int, F_int, mobility_mp, k_mp, reduction_factor, k_L
	Variable i
	Variable index = 0
	for(i=0;i<numpnts(illuminations);i+=1)
		if(rejection_mask[i]==1)
			continue
		endif
		// Gather JV information
		SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV:$(illuminations[i])
		Wave J_gen
		Wave voltage_effective
		Wave J_br
		NVAR J_sc
		NVAR V_oc
		NVAR V_0
		NVAR FF
		NVAR V_mp
		NVAR J_mp
		NVAR V_oc
		NVAR P_max
		NVAR I_suns_val = $("I_suns")
		SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance:$("CV_"+illuminations[i])
		NVAR n_mp
		Wave voltage_applied
		Wave voltage_cor
		Wave V_ints = $("V_int")
		Wave F_ints = $("F_int")
		Wave k_br
		Wave mu_eff
		// Calculate device metrics
		I_suns[index] = {I_suns_val}
		n_mps[index] = {n_mp}
		V_int[index] = {interp(V_mp,voltage_applied,V_ints)}
		F_int[index] = {V_int[index]/Active_thickness_cm}
		G_sc[index] = {1e-3*J_gen(0)/(Active_thickness_cm*1.60217662e-19)}
		G_mp[index] = {1e-3*J_gen(V_mp)/(Active_thickness_cm*1.60217662e-19)}
		mobility_mp[index] = {interp(V_mp,voltage_applied,mu_eff)}
		k_mp[index] = {interp(V_mp,voltage_applied,k_br)}
		J_scs[index] = {J_sc}
		FFs[index] = {FF}
		V_ocs[index] = {V_oc}
		PCEs[index] = {P_max/I_suns_val}
		tau_ex[index] = {Active_thickness_cm^2/(2*mobility_mp[index]*V_int[index])}
		tau_rec[index] = {1/(k_mp[index]*n_mp)}
		theta[index] = {(k_mp[index]*G_mp[index]*Active_thickness_cm^4)/(V_int[index]^2*mobility_mp[index]^2)}
		theta_alpha[index] = {sqrt((k_mp[index]*G_mp[index]*Active_thickness_cm^4)/(4*mobility_mp[index]^2*(8.6173e-5*300)^2))}
		gamma_dir[index] = {(k_mp[index]^0.8*(-J_sc)^0.5*Active_thickness_cm^3.5)/(V_oc^2*mobility_mp[index]^2)}
		k_L[index] = {(2*1.602177e-19*mobility_mp[index])/(8.854187e-12*1e-2*epsilon)}
		reduction_factor[index] = {k_mp[index]/k_L[index]}
		index += 1
	endfor
	// Construct J_gen and J_br matrices
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):IPDA_Results
	Make/N=(numpnts(J_br),numpnts(illuminations))/O J_br_intensity, J_gen_intensity
	index = 0
	for(i=0;i<numpnts(illuminations);i+=1)
		if(rejection_mask[i]==1)
			continue
		endif
		SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV:$(illuminations[i])
		Wave J_gen
		Wave J_br
		J_gen_intensity[][index] = J_gen[p]
		J_br_intensity[][index] = J_br[p]
		index += 1
	endfor
	SetScale/P x -4,DimDelta(J_br,0),"", J_gen_intensity
	SetScale/P x -4,DimDelta(J_br,0),"", J_br_intensity
	// Calculate current and FF without field dependent generation and without series resistance 
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):IPDA_Results	
	Wave I_suns
	Make/D/O/N=0 FF_brs, FF_0s
	index = 0
	for(i=0;i<numpnts(illuminations);i+=1)
		if(rejection_mask[i]==1)
			continue
		endif
		SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV:dark
		NVAR R_s
		SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV:$(illuminations[i])
		NVAR V_0
		Wave J_gen
		Duplicate/O J_gen, $"J_br_nfg", $"J_nfg"
		Wave J_br_nfg
		Wave J_nfg
		Variable J_gen_nfg = J_gen(V_0-2)
		Variable j
		for(j=0;j<numpnts(J_gen);j+=1)
			Make/N=(numpnts(I_suns))/O gen_wave, br_wave
			gen_wave = J_gen_intensity[j][p]
			br_wave = J_br_intensity[j][p]
			Variable intensity_target = interp(J_gen_nfg,gen_wave,I_suns)
			J_br_nfg[j] = interp(intensity_target,I_suns,br_wave)
		endfor
		J_nfg = J_br_nfg - J_gen_nfg
		// Calculate JV metrics
		FEDMS_CalculatePVJVResults(J_nfg)
		Wave W_result
		FF_brs[index] = {W_result[5]}
		Variable r_s_norm = R_s*-W_result[0]*Device_area_cm2/(1000*W_result[1])
		FF_0s[index] = {(FF_brs[index]-(r_s_norm^2/5.4))/(1-1.1*r_s_norm)}
		index += 1
		KillWaves/Z gen_wave br_wave W_result fit_J_nfg W_fitConstants
	endfor
	// Create results summary table
	if(hide_table!=1)
		SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):IPDA_Results	
		Edit/W=(400,400,1480,630) I_suns,G_mp,J_scs,V_ocs,FFs,PCEs,n_mps,k_mp,mobility_mp,reduction_factor,tau_rec,tau_ex,theta as device_name+"_IPDA_Results"
		ModifyTable width=75,format(Point)=1 
	endif
	SetDataFolder $(original_data_folder)
End

Function FEDMS_CalculateIPDAResultsGUI()
	String sample_name = FEDMS_ChooseDevice("Impedance")
	if(StringMatch(sample_name,""))
		return NaN
	endif
	Print "•FEDMS_CalculateIPDAResults(\""+sample_name+"\")"
	FEDMS_CalculateIPDAResults(sample_name)
End

Function FEDMS_CalculatePhotocurrent(device_name,measurement_name)
	String device_name
	String measurement_name
	String original_data_folder = GetDataFolder(1)
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[strlen(device_name)-1]
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV
	// Calculate photocurrent
	// Get appropriate dark current data
	if(StringMatch(measurement_name,"*LED*"))
		FEDMS_AnalyzeJV(device_name,"JV","LEDdark")
		SetDataFolder :LEDdark
	else
		FEDMS_AnalyzeJV(device_name,"JV","dark")
		SetDataFolder :dark
	endif
	Wave J_dark = $("J_avg")
	// Get the light current data
	FEDMS_AnalyzeJV(device_name,"JV",measurement_name)
	SetDataFolder ::$measurement_name
	Wave J_light = $("J_avg")
	Duplicate/O J_light photocurrent
	photocurrent = -(J_light-J_dark)
	// Calculate V_0 and the effective voltage waveform
	WaveStats/Q photocurrent
	FindLevel/Q/R=(0,pnt2x(photocurrent,V_npnts-1)) photocurrent,0
	Reverse photocurrent
	Variable/G/D V_0 = V_LevelX
	if(numtype(V_0)==2)
		WaveStats/Q photocurrent
		V_0 = V_minloc
	endif
	Duplicate/O photocurrent voltage_effective
	voltage_effective = V_0 - x
	SetDataFolder original_data_folder
End

Function FEDMS_CalculatePVJVResults(J_wave)
	Wave J_wave
	// Calculate the open-circuit voltage and short-circuit current
	FindLevel/Q J_wave,0
	Variable V_oc = V_LevelX
	Variable J_sc = J_wave(0)
	// Calculate the maximum power point data
	//  Determine simple first guess
	Duplicate/O J_wave power
	power = -J_wave*x
	WaveStats/Q power
	Variable P_max = V_max
	Variable V_mp = V_maxloc
	//  Fit current curve near estimated max power to obtain more accurate value
	if(V_oc > 0.2)
		Variable index_start = x2pnt(J_wave,V_maxloc-0.1)
		Variable index_end = x2pnt(J_wave,V_maxloc+0.05)
		Variable offset = pnt2x(J_wave,index_start)
		K0 = J_wave(pnt2x(J_wave,index_start));K1 = 0.1;K2 = -0.05;
		CurveFit/G/N/NTHR=0/W=2/Q/K={offset} exp_XOffset  J_wave[index_start,index_end] /D 
		Wave W_coef
		Duplicate/O J_wave deriv
		deriv = W_coef[0]-(W_coef[1]/W_coef[2])*x*exp(-(x-offset)/W_coef[2])+W_coef[1]*exp(-(x-offset)/W_coef[2])
		FindLevel/Q deriv,0
		V_mp = V_LevelX
		P_max = -W_coef[0]*V_mp-W_coef[1]*V_mp*exp(-(V_mp-offset)/W_coef[2])
	endif
	Variable J_mp = J_wave(V_mp)
	Variable FF = J_mp*V_mp/(J_sc*V_oc)
	Make/N=6/D/O W_result
	W_result[0] = J_sc
	W_result[1] = V_oc
	W_result[2] = P_max
	W_result[3] = J_mp
	W_result[4] = V_mp
	W_result[5] = FF
	KillWaves/Z deriv, power, $("fit_J_wave") W_coef $"W_sigma"
End

Function/S FEDMS_ChooseDevice(keyword)
	String keyword
	String original_folder = GetDataFolder(1)
	// Build the sample list
	SetDataFolder root:FlexEDMS
	String device_list = ""
	DFREF dfr1 = GetDataFolderDFR()
	Variable N_samples = CountObjectsDFR(dfr1,4)
	String substrate_name
	String device_name
	Variable i
	Variable j
	for(i=0;i<N_samples;i+=1)
		substrate_name = GetIndexedObjNameDFR(dfr1,4,i)
		if(StringMatch("Analysis",substrate_name)==0)
			SetDataFolder root:FlexEDMS:$(substrate_name)
			DFREF dfr2 = GetDataFolderDFR()
			Variable N_devices = CountObjectsDFR(dfr2,4)
			for(j=0;j<N_devices;j+=1)
				device_name = GetIndexedObjNameDFR(dfr2,4,j)
				SetDataFolder root:FlexEDMS:$(substrate_name):$(device_name)
				if(DataFolderExists(":"+keyword))
					device_list = AddListItem(substrate_name+device_name,device_list)
				endif
			endfor
		endif
	endfor
	String sample_name
	// Prompt user to choose the sample
	Prompt sample_name, "Choose the sample name:", popup, device_list
	DoPrompt "Make Selection",sample_name
	// User cancelled operation
	if(V_flag==1)
		SetDataFolder original_folder
		return ""
	endif
	SetDataFolder original_folder
	return sample_name
End

Function/S FEDMS_ChooseMeasurement(device_name,measurement_type)
	String device_name
	String measurement_type
	String original_folder = GetDataFolder(1)
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[strlen(device_name)-1]
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):$(measurement_type)
	String measurement_name
	String measurement_list = ""
	DFREF dfr1 = GetDataFolderDFR()
	Variable i
	Variable N_measurements = CountObjectsDFR(dfr1,4)
	for(i=0;i<N_measurements;i+=1)
		measurement_name = GetIndexedObjNameDFR(dfr1,4,i)
		measurement_list = AddListItem(measurement_name,measurement_list)
	endfor
	// Prompt user to choose the measurement
	Prompt measurement_name, "Choose the measurement condition:", popup, measurement_list
	DoPrompt "Make Selection",measurement_name
	// User cancelled operation
	if(V_flag==1)
		SetDataFolder original_folder
		return ""
	endif
	SetDataFolder original_folder
	return measurement_name
End

Function/S FEDMS_ChooseTOFCarrierType(sample_name)
	String sample_name
	String carrier_type
	String original_folder = GetDataFolder(1)
	// Build the carrier type list
	SetDataFolder root:FlexEDMS:$(sample_name):$("Time of Flight")
	String carrier_type_list = ""
	DFREF dfr2 = GetDataFolderDFR()
	Variable N_types = CountObjectsDFR(dfr2,4)
	Variable i
	String folder_name
	for(i=0;i<N_types;i+=1)
		folder_name = GetIndexedObjNameDFR(dfr2,4,i)
		carrier_type_list = AddListItem(folder_name,carrier_type_list)
	endfor
	Prompt carrier_type, "Choose the carrier type:", popup, carrier_type_list
	DoPrompt "Make Selection",carrier_type
	if(V_flag==1)
		SetDataFolder original_folder
		return ""
	endif
	SetDataFolder original_folder
	return carrier_type
End

Function/S FEDMS_ChooseSample()
	String original_folder = GetDataFolder(1)
	// Build the sample list
	String sample_name
	SetDataFolder root:FlexEDMS
	String sample_list = ""
	DFREF dfr1 = GetDataFolderDFR()
	Variable N_samples = CountObjectsDFR(dfr1,4)
	String folder_name
	Variable i
	for(i=0;i<N_samples;i+=1)
		folder_name = GetIndexedObjNameDFR(dfr1,4,i)
		sample_list = AddListItem(folder_name,sample_list)
	endfor
	// Prompt user to choose the sample
	Prompt sample_name, "Choose the sample name:", popup, sample_list
	DoPrompt "Make Selections",sample_name
	// User cancelled operation
	if(V_flag==1)
		SetDataFolder original_folder
		return ""
	endif
	SetDataFolder original_folder
	return sample_name
End

Function/S FEDMS_ChooseSubfolder(folder_path)
	String folder_path
	String original_folder = GetDataFolder(1)
	// Build the subfolder list
	SetDataFolder folder_path
	String subfolder_list = ""
	DFREF dfr1 = GetDataFolderDFR()
	Variable N_folders = CountObjectsDFR(dfr1,4)
	String subfolder_name
	Variable i
	for(i=0;i<N_folders;i+=1)
		subfolder_name = GetIndexedObjNameDFR(dfr1,4,i)
		subfolder_list = AddListItem(subfolder_name,subfolder_list)
	endfor
	// Prompt user to choose the subfolder
	Prompt subfolder_name, "Choose the subfolder:", popup, subfolder_list
	DoPrompt "Make Selection",subfolder_name
	// User cancelled operation
	if(V_flag==1)
		SetDataFolder original_folder
		return ""
	endif
	SetDataFolder original_folder
	return subfolder_name
End

Function/S FEDMS_ChooseTOFSample()
	String original_folder = GetDataFolder(1)
	// Build the sample list
	String sample_name
	SetDataFolder root:FlexEDMS
	String sample_list = ""
	DFREF dfr1 = GetDataFolderDFR()
	Variable N_samples = CountObjectsDFR(dfr1,4)
	String folder_name
	Variable i
	for(i=0;i<N_samples;i+=1)
		folder_name = GetIndexedObjNameDFR(dfr1,4,i)
		SetDataFolder :$(folder_name)
		if(DataFolderExists("Time of Flight"))
			sample_list = AddListItem(folder_name,sample_list)
		endif
		SetDataFolder ::
	endfor
	// Prompt user to choose the sample
	Prompt sample_name, "Choose the sample name:", popup, sample_list
	DoPrompt "Make Selections",sample_name
	// User cancelled operation
	if(V_flag==1)
		SetDataFolder original_folder
		return ""
	endif
	SetDataFolder original_folder
	return sample_name
End

Function/S FEDMS_DoOpenMultiFileDialog()
	Variable refNum
	String message = "Select one or more files"
	String outputPaths = ""
	String fileFilters = "Data Files (*.txt,*.dat,*.csv):.txt,.dat,.csv;"
	fileFilters += "All Files:.*;"
	Open /D /R /MULT=1 /F=fileFilters /M=message refNum
	outputPaths = S_fileName
	if (strlen(outputPaths) == 0)
		return ""
	endif
	return outputPaths		// Will be empty if user canceled
End

Function FEDMS_EnterSampleInfo()
	String original_folder = GetDataFolder(1)
	// Build the sample list
	SetDataFolder root:FlexEDMS
	String sample_list = ""
	DFREF dfr1 = GetDataFolderDFR()
	Variable N_samples = CountObjectsDFR(dfr1,4)
	String folder_name
	Variable i
	for(i=0;i<N_samples;i+=1)
		folder_name = GetIndexedObjNameDFR(dfr1,4,i)
		if(StringMatch(folder_name,"Analysis")==0)
			sample_list = AddListItem(folder_name,sample_list)
		endif
	endfor
	sample_list = AddListItem("New Sample",sample_list)
	// Prompt user to create new or choose existing sample
	String sample_name
	Prompt sample_name, "Create new sample or choosing existing sample:", popup, sample_list
	DoPrompt "Make Selection",sample_name
	// User cancelled operation
	if(V_flag==1)
		SetDataFolder original_folder
		return NaN
	endif
	// Build sample info form
	String fab_date
	String fab_persons
	String sample_comp
	String ratio
	String solvent
	Variable concentration
	Variable speed
	Variable sample_thickness
	Variable annealing_temp
	Variable annealing_time
	Variable N_devices
	Variable device_area
	String comments
	Prompt sample_name, "Enter the sample name:"
	Prompt fab_date, "Enter the fabrication date (mm/dd/yyyy):"
	Prompt fab_persons, "Enter the fabrication person(s):"
	Prompt sample_comp, "Enter the composition of the semiconductor layer:"
	Prompt ratio, "Enter the blend ratio (wt:wt):"
	Prompt solvent, "Enter the name of the casting solvent:"
	Prompt concentration, "Enter the concentration of the solution (mg/mL):"
	Prompt speed, "Enter the coating speed (mm/s or rpm):"
	Prompt sample_thickness, "Enter the semiconductor layer thickness (cm):"
	Prompt annealing_temp, "Enter the annealing temperature (C):"
	Prompt annealing_time, "Enter the annealing time (min):"
	Prompt N_devices, "Enter the number of devices on the substrate:"
	Prompt device_area, "Enter the active area of the devices (cm^2):"
	Prompt comments, "Enter any additional comments and notes:"
	// Load existing sample data into form variables if they exist
	if(StringMatch("New Sample",sample_name)==0)
		SetDataFolder :$(sample_name)
		SVAR Fabrication_date
		SVAR Fabrication_persons
		SVAR Sample_composition
		SVAR Blend_ratio
		SVAR Casting_solvent
		NVAR Annealing_temperature_C
		NVAR Annealing_time_min
		NVAR Active_thickness_cm
		NVAR Device_area_cm2
		SVAR Sample_comments
		if(Exists("Fabrication_date"))
			fab_date = Fabrication_date
		endif
		if(Exists("Fabrication_persons"))
			fab_persons = Fabrication_persons
		endif
		if(Exists("Sample_composition"))
			sample_comp = Sample_composition
		endif
		if(Exists("Blend_ratio"))
			ratio = Blend_ratio
		endif
		if(Exists("Casting_solvent"))
			
		endif
		if(Exists("Active_thickness_cm"))
			sample_thickness = Active_thickness_cm
		endif
		if(Exists("Device_area_cm2"))
			device_area = Device_area_cm2
		endif
		if(Exists("Sample_comments"))
			comments = Sample_comments
		endif
		if(Exists("Annealing_temperature_C"))
			annealing_temp = Annealing_temperature_C
		endif
		if(Exists("Annealing_time_min"))
			annealing_time = Annealing_time_min
		endif
	else
		sample_name = ""
	endif
	// Prompt user to enter/update sample info
	Variable input_needed = 1
	do
		// First prompt
		DoPrompt "Enter Sample Info", sample_name, fab_date, fab_persons, sample_comp, ratio, solvent, concentration
		// Check for user cancel
		if(V_flag==1)
			break
		endif
		// Second prompt
		DoPrompt "Enter Sample Info", speed, sample_thickness, annealing_temp, annealing_time, N_devices, device_area, comments
		// Check for user cancel
		if(V_flag==1)
			break
		endif
		// Check if required parameters have been completed
		if(!sample_thickness>0 || StringMatch(sample_name,""))
			DoAlert 0, "You must at least specify a sample name and a layer thickness value."
			continue
		// Check for sample name conflict
		elseif(FindListItem(sample_name,sample_list)!=-1)
			DoAlert 0, "A sample with that name already exists in your database, enter a different name."
		else
			input_needed = 0
		endif
	while(input_needed==1)
	// User cancelled operation
	if(V_flag==1)
		SetDataFolder original_folder
		return NaN
	endif
	// Save user input into sample folder
	if(!DataFolderExists(sample_name))
		NewDataFolder $(sample_name)
	endif
	SetDataFolder :$(sample_name)
	String/G Fabrication_date = fab_date
	String/G Fabrication_persons = fab_persons
	String/G Sample_composition = sample_comp
	String/G Blend_ratio = ratio
	String/G Casting_solvent = solvent
	Variable/G Solution_concentration = concentration
	Variable/G Coating_speed = speed
	Variable/G Annealing_temperature_C = annealing_temp
	Variable/G Annealing_time_min = annealing_time
	Variable/G Active_thickness_cm = sample_thickness
	Variable/G Device_area_cm2 = device_area
	String/G Sample_comments = comments
	if(N_devices>0)
		for(i=0;i<N_devices;i+=1)
			NewDataFolder $(num2char(i+65))
		endfor
	endif
	SetDataFolder original_folder
End

Function FEDMS_FilterMobilityData(sample_name,carrier_type,dispersion_limit)
	String sample_name
	String carrier_type
	Variable dispersion_limit
	String original_folder = GetDataFolder(1)
	SetDataFolder root:FlexEDMS:$(sample_name):$("Time of Flight"):$(carrier_type)
	Wave dispersion
	Wave temperature_K
	Wave mobility_geo
	Wave rejection_mask
	Make/N=1/D/O mobility_filtered, temp_filtered, dispersion_filtered
	Variable i
	Variable count = 0
	for(i=0;i<numpnts(mobility_geo);i+=1)
		if(dispersion[i]<dispersion_limit && rejection_mask[i]==0)
			dispersion_filtered[count] = {dispersion[i]}
			temp_filtered[count] = {temperature_K[i]}
			mobility_filtered[count] = {mobility_geo[i]}
			count += 1
		endif
	endfor
	SetDataFolder original_folder
End

Function FEDMS_FitInductance(w,f) : FitFunc
	Wave w
	Variable f
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(f) = 2*PI*f*L-(1/(2*PI*f*C))
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ f
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w[0] = L
	//CurveFitDialog/ w[1] = C
	return 2*PI*f*w[0]-(1/(2*PI*f*w[1]))
End

Function FEDMS_FitMobilityField(w,F_sqrt) : FitFunc
	Wave w
	Variable F_sqrt
	//CurveFitDialog/ These comments were created by the Curve Fitting dialog. Altering them will
	//CurveFitDialog/ make the function less convenient to work with in the Curve Fitting dialog.
	//CurveFitDialog/ Equation:
	//CurveFitDialog/ f(F_sqrt) = mu_0*exp(beta*F_sqrt)
	//CurveFitDialog/ End of Equation
	//CurveFitDialog/ Independent Variables 1
	//CurveFitDialog/ F_sqrt
	//CurveFitDialog/ Coefficients 2
	//CurveFitDialog/ w[0] = mu_0
	//CurveFitDialog/ w[1] = beta
	return w[0]*exp(w[1]*F_sqrt)
End

Function FEDMS_LoadImpedanceFile(sample_name,fullpathname,file_format_in,device_area,active_thickness)
	String sample_name
	String fullpathname
	String file_format_in
	Variable device_area
	Variable active_thickness
	// Parse filename from full pathname
	String filename = StringFromList(ItemsInList(fullpathname,":")-1,fullpathname,":")
	// Parse filename to extract device name
	String file_device_name = StringFromList(1,filename,"_")
	// Check that filename and sample name match
	String file_batch_name = StringFromList(0,filename,"_")
	String file_sample_name = file_batch_name+"_"+file_device_name
	file_sample_name = RemoveEnding(file_sample_name)
	if(!StringMatch(sample_name,file_sample_name))
		Print "Error! The file sample name does not match the sample you selected."
		return NaN
	endif
	String original_folder = GetDataFolder(1)
	NewDataFolder/S/O root:FlexEDMS
	NewDataFolder/S/O $(sample_name)
	String/G File_format = file_format_in
	Variable/G Device_area_cm2 = device_area
	Variable/G Active_thickness_cm = active_thickness
	// Parse filename to extract measurement conditions
	String measurement_type = StringFromList(2,filename,"_")
	String measurement_name = StringFromList(3,filename,"_")
	measurement_name = RemoveEnding(measurement_name,".txt")
	// Create Subfolders
	String device_letter = file_device_name[strlen(file_device_name)-1]
	NewDataFolder/O/S $device_letter
	NewDataFolder/O/S $"Impedance"
	// If not a dark measurement, add to illuminations wave
	if(!StringMatch(measurement_name,"*dark*"))
		Wave/T illuminations = $("illuminations")
		// Create illuminations wave if it does not exists
		if(!WaveExists(illuminations))
			Make/T/N=1 illuminations
			illuminations[0] = measurement_name
		// Check if illumination value already exists
		else
			Variable item_exists = 0
			Variable i
			for(i=0;i<numpnts(illuminations);i+=1)
				if(StringMatch(illuminations[i],measurement_name))
					item_exists = 1
				endif
			endfor
			// If it does not exist, expand the wave and add it
			if(item_exists==0)
				illuminations[numpnts(illuminations)] = {measurement_name}
			endif
		endif
	endif
	if(StringMatch(measurement_type,"*Cf*"))
		NewDataFolder/O/S $("Cf_"+measurement_name)
	elseif(StringMatch(measurement_type,"*CV*"))
	NewDataFolder/O/S $("CV_"+measurement_name)
	endif
	if(StringMatch(file_format,"Richter Lab"))
		// Load file into text wave
		LoadWave/A/J/Q/K=2/V={""," $",0,0} fullpathname
		Wave/T wave0
		// Parse file header for measurement date and time
		String date_str = StringFromList(2,wave0[0]," ")
		date_str = RemoveEnding(date_str)
		String day = StringFromList(0,date_str,".")
		String month = StringFromList(1,date_str,".")
		String year = StringFromList(2,date_str,".")
		String/G Measurement_date = month+"/"+day+"/"+year
		String/G Measurement_time = (StringFromList(3,wave0[0]," "))
		Variable/G Amplitude_V = str2num(StringFromList(1,wave0[1],"="))
		// Load Data Waves
		LoadWave/A/J/Q/O/K=0/L={2,3,0,0,0} fullpathname
		// Load C-f scan data
		if(StringMatch(measurement_type,"*Cf*"))
			// Split into multiple waves based on the dc bias
			// Create list of biases tested
			Wave wave2
			Make/N=0/O voltage_list
			Make/N=0/O start_index_list
			Make/N=0/O end_index_list
			Variable current_val = wave2[0]
			voltage_list[0] = {current_val}
			start_index_list[0] = {0}
			for(i=1;i<numpnts(wave2);i+=1)
				if(wave2[i]!=current_val)
					current_val = wave2[i]
					end_index_list[numpnts(end_index_list)] = {i-1}
					start_index_list[numpnts(start_index_list)] = {i}
					voltage_list[numpnts(voltage_list)] = {current_val}
				endif
			endfor
			end_index_list[numpnts(end_index_list)] = {numpnts(wave2)-1}
			// Split the temporary waves
			Duplicate/O/R=[start_index_list[0],end_index_list[0]] $("wave1") frequency
			Reverse frequency
			SetScale/I x 0,numpnts(frequency)-1,"", frequency
			for(i=0;i<numpnts(voltage_list);i+=1)
				Duplicate/O/R=[start_index_list[i],end_index_list[i]] $("wave3") $("Z_real_"+num2str(voltage_list[i])+"V")
				Duplicate/O/R=[start_index_list[i],end_index_list[i]] $("wave4") $("Z_imag_"+num2str(voltage_list[i])+"V")
				Reverse $("Z_real_"+num2str(voltage_list[i])+"V")
				SetScale/I x 0,numpnts(frequency)-1,"", $("Z_real_"+num2str(voltage_list[i])+"V")
				Reverse $("Z_imag_"+num2str(voltage_list[i])+"V")
				SetScale/I x 0,numpnts(frequency)-1,"", $("Z_imag_"+num2str(voltage_list[i])+"V")
			endfor
		// Load C-V scan data
		elseif(StringMatch(measurement_type,"*CV*"))
			Variable/G Freq = str2num(StringFromList(0,StringFromList(1,wave0[1],"=")," "))
			Duplicate/O $("wave1") voltage_applied
			Duplicate/O $("wave2") Z_real
			Duplicate/O $("wave3") Z_imag
		endif
	elseif(StringMatch(file_format,"Hersam Lab"))
		// Load file into text wave
		LoadWave/A/J/Q/K=2/V={""," $",0,0} fullpathname
		Wave/T wave0
		// Parse file header for measurement date and time
		date_str = StringFromList(ItemsInList(wave0[3]," ")-1,wave0[3]," ")
		date_str = ReplaceString("-",date_str,"/")
		String/G Measurement_date = date_str
		String/G Measurement_time = StringFromList(ItemsInList(wave0[4]," ")-1,wave0[4]," ")
		// Find end of header section
		FindValue/TEXT="End Comments" wave0
		Variable data_start_index = V_value+1
		// Load Data Waves
		LoadWave/A/J/D/Q/O/K=0/L={data_start_index-2,data_start_index,0,0,0} fullpathname
		Wave wave2
		Variable/G Amplitude_V = wave2[0]
		// Load dark C-f scan data
		if(StringMatch(measurement_type,"*Cf*"))
			// Split into multiple waves based on the dc bias
			// Create list of biases tested
			Wave wave2, wave3
			Make/N=0/O voltage_list
			voltage_list[0] = {wave3[0]}
			Duplicate/O $("wave1") frequency
			Duplicate/O $("wave5") $("Z_real_"+num2str(voltage_list[0])+"V")
			Duplicate/O $("wave6") $("Z_imag_"+num2str(voltage_list[0])+"V")
		// Load C-V scan data
		elseif(StringMatch(measurement_type,"*CV*"))
			Wave wave1
			Variable/G Freq = wave1[0]
			Duplicate/O $("wave3") voltage_applied
			Duplicate/O $("wave5") Z_real
			Duplicate/O $("wave6") Z_imag
		endif
	endif
	// Clean up
	KillWaves/Z $"wave0", $"wave1", $"wave2", $"wave3", $"wave4", $"wave5", $"wave6", $"wave7", $"wave8", $"wave9"
	KillWaves/Z voltage_list, start_index_list, end_index_list
	SetDataFolder original_folder
End

Function FEDMS_LoadImpedanceFolder(path_str, measurement_persons, comments, file_format, device_area, active_thickness)
	String path_str
	String measurement_persons
	String comments
	String file_format
	Variable device_area
	Variable active_thickness
	String original_folder = GetDataFolder(1)
	// Get list of all txt files in the folder
	NewPath/O/Q folder_path, path_str
	String file_list
	if(StringMatch(file_format,"Richter Lab"))
		file_list = IndexedFile(folder_path,-1,".txt")
	elseif(StringMatch(file_format,"Hersam Lab"))
		file_list = IndexedFile(folder_path,-1,".z")
	endif
	// Filter list for impedance measurement data
	Variable i
	for(i=0;i<ItemsInList(file_list,";");i+=1)
		String file_name = StringFromList(i,file_list,";")
		// Filter out non-impedance files
		if(!StringMatch(file_name,"*impedance*"))
			file_list = RemoveListItem(i,file_list)
			i -= 1
		endif
	endfor
	file_list = SortList(file_list,";",16)
	Variable file_count = ItemsInList(file_list)
	// Load data from each file
	Print "Loading "+num2str(ItemsInList(file_list,";"))+" impedance data files..."
	for(i=0;i<ItemsInList(file_list,";");i+=1)
		file_name = StringFromList(i,file_list,";")
		String fullpathname = path_str+file_name
		String sample_name = StringFromList(0,file_name,"_")+"_"+StringFromList(1,file_name,"_")
		sample_name = RemoveEnding(sample_name)
		FEDMS_LoadImpedanceFile(sample_name,fullpathname,file_format,device_area,active_thickness)
	endfor
	SetDataFolder original_folder
End

Function FEDMS_LoadImpedanceFolderGUI()
	String measurement_persons
	String measurement_comments
	String file_format
	Variable device_area
	Variable active_thickness
	Prompt measurement_persons, "Enter the measurement person(s):"
	Prompt measurement_comments, "Enter any additional measurement comments:"
	String format_list = "Richter Lab;Hersam Lab"
	Prompt file_format, "Choose the file format:", popup, format_list
	Prompt device_area, "Enter the device area (cm^2):"
	Prompt active_thickness, "Enter the active layer thickness (cm):"
	DoPrompt "Enter Measurement Info", measurement_persons, measurement_comments, file_format, device_area, active_thickness
	if(V_flag!=0)
		return NaN
	endif
	NewPath/O/Q folder_path
	if(V_flag!=0)
		return NaN
	endif
	if(device_area<=0)
		Print "Error! invalid device area entered."
		return NaN
	endif
	if(device_area<=0)
		Print "Error! invalid active layer thickness entered."
		return NaN
	endif
	PathInfo folder_path
	Print "•FEDMS_LoadImpedanceFolder(\""+S_path+"\",\""+measurement_persons+"\",\""+measurement_comments+"\",\""+file_format+"\","+num2str(device_area)+","+num2str(active_thickness)+")"
	FEDMS_LoadImpedanceFolder(S_path, measurement_persons, measurement_comments, file_format,device_area,active_thickness)
End

Function FEDMS_LoadJVFile(sample_name,fullpathname,persons,comments,file_format_in,device_area,mismatch)
	String sample_name
	String fullpathname
	String persons
	String comments
	String file_format_in
	Variable device_area
	Variable mismatch
	// Parse filename from full pathname
	String filename = StringFromList(ItemsInList(fullpathname,":")-1,fullpathname,":")
	// Parse filename to extract device name
	String file_device_name = StringFromList(1,filename,"_")
	// Check that filename and sample name match
	String file_batch_name = StringFromList(0,filename,"_")
	String file_sample_name = file_batch_name+"_"+file_device_name
	file_sample_name = RemoveEnding(file_sample_name)
	if(!StringMatch(sample_name,file_sample_name))
		Print "Error! The file sample name does not match the sample you selected."
		return NaN
	endif
	String original_folder = GetDataFolder(1)
	NewDataFolder/O/S root:FlexEDMS
	NewDataFolder/O/S $(sample_name)
	String/G File_format = file_format_in
	Variable/G Device_area_cm2 = device_area
	// Parse filename to extract measurement conditions
	String measurement_type = RemoveEnding(StringFromList(2,filename,"_"),".txt")
	String measurement_name = RemoveEnding(StringFromList(3,filename,"_"),".txt")
	// Create Subfolders
	String device_letter = file_device_name[strlen(file_device_name)-1]
	NewDataFolder/O/S $device_letter
	NewDataFolder/O/S $measurement_type
	// Enter measurement info
	String/G Measurement_persons = persons
	String/G Measurement_comments = comments
	Variable/G Mismatch_factor = mismatch
	// If not a dark measurement, add to illuminations wave
	if(!StringMatch(measurement_name,"*dark*"))
		Wave/T illuminations
		// Create illuminations wave if it does not exist
		if(!WaveExists(illuminations))
			Make/T/N=1 illuminations
			illuminations[0] = measurement_name
		// Check if illumination value already exists
		else
			Variable item_exists = 0
			Variable i
			for(i=0;i<numpnts(illuminations);i+=1)
				if(StringMatch(illuminations[i],measurement_name))
					item_exists = 1
				endif
			endfor
			// If it does not exist, expand the wave and add it
			if(item_exists==0)
				illuminations[numpnts(illuminations)] = {measurement_name}
			endif
		endif
	endif
	NewDataFolder/O/S $measurement_name
	if(StringMatch(file_format,"Richter Lab"))
		// Load file into text wave
		LoadWave/A/J/Q/K=2/V={""," $",0,0} fullpathname
		Wave/T wave0
		// Parse file header for measurement date and time
		String/G Measurement_date = StringFromList(1,StringFromList(1,wave0[0],";"),":")
		String am_pm = StringFromList(4,wave0[0]," ")
		am_pm = RemoveEnding(am_pm)
		String time_str = StringFromList(3,wave0[0]," ")
		if(StringMatch(am_pm,"PM"))
			Variable hour = str2num(StringFromList(0,time_str,":"))
			String minute_str = StringFromList(1,time_str,":")
			hour += 12
			time_str = num2str(hour)+":"+minute_str
		endif
		String/G Measurement_time = time_str
		// Load Data Waves
		LoadWave/A/J/D/Q/O/K=0/O/L={1,2,0,0,0} fullpathname
		Duplicate/O $("wave1") voltage
		Duplicate/O $("wave2") current
	elseif(StringMatch(file_format,"Hersam Lab"))
		// Load Data Waves
		LoadWave/A/J/D/Q/O/K=0/O/L={0,1,0,0,0} fullpathname
		Duplicate/O $("wave1") voltage
		Duplicate/O $("wave0") current
		// round voltage values to nearest mV
		for(i=0;i<numpnts(voltage);i+=1)
			voltage[i] = round(voltage[i]*1000)/1000
		endfor
	endif
	// Check for inverted measurement
	Duplicate voltage voltage_abs
	voltage_abs = abs(voltage)
	WaveStats/Q voltage_abs
	if(current[V_maxloc]>0)
		voltage *= -1
		current *= -1
	endif
	// Clean up
	KillWaves/Z $"wave0" $"wave1" $"wave2" $"wave3" $"wave4" voltage_abs
	SetDataFolder original_folder
End

Function FEDMS_LoadJVFolder(path_str, measurement_persons, comments, file_format, device_area,mismatch_factor)
	String path_str
	String measurement_persons
	String comments
	String file_format
	Variable device_area // cm^2
	Variable mismatch_factor
	String original_folder = GetDataFolder(1)
	// Get list of all txt files in the folder
	NewPath/O/Q folder_path , path_str
	String file_list = IndexedFile(folder_path,-1,".txt")
	// Filter list for J-V measurement data
	Variable i
	for(i=0;i<ItemsInList(file_list,";");i+=1)
		String file_name = StringFromList(i,file_list,";")
		// Filter out non-JV files
		if(!StringMatch(file_name,"*JV*"))
			file_list = RemoveListItem(i,file_list)
			i -= 1
			continue
		endif
		// Filter out summary files
		if(StringMatch(file_name,"*Summary*"))
			file_list = RemoveListItem(i,file_list)
			i -= 1
			continue
		endif
	endfor
	file_list = SortList(file_list,";",16)
	Variable file_count = ItemsInList(file_list)
	// Load data from each file
	Print "Loading "+num2str(ItemsInList(file_list,";"))+" JV data files..."
	for(i=0;i<ItemsInList(file_list,";");i+=1)
		file_name = StringFromList(i,file_list,";")
		String fullpathname = path_str+file_name
		String sample_name = StringFromList(0,file_name,"_")+"_"+StringFromList(1,file_name,"_")
		sample_name = RemoveEnding(sample_name)
		FEDMS_LoadJVFile(sample_name,fullpathname,measurement_persons,comments,file_format,device_area,mismatch_factor)
	endfor
	SetDataFolder original_folder
End

Function FEDMS_LoadJVFolderGUI()
	String original_folder = GetDataFolder(1)
	String measurement_persons
	String measurement_comments
	String file_format
	Variable device_area = 0
	Variable mismatch_factor = 1
	Prompt measurement_persons, "Enter the measurement person(s):"
	Prompt measurement_comments, "Enter any additional measurement comments:"
	String format_list = "Richter Lab;Hersam Lab"
	Prompt file_format, "Choose the file format:", popup, format_list
	Prompt device_area, "Enter the device area (cm^-2):"
	Prompt mismatch_factor, "Enter the 1 sun mismatch factor:"
	DoPrompt "Enter Measurement Info", measurement_persons, measurement_comments, file_format, device_area, mismatch_factor
	if(V_flag==1)
		return NaN
	endif
	if(device_area<=0)
		Print "Error! Invalid device area."
		return NaN
	endif
	NewPath/O/Q folder_path
	if(V_flag!=0)
		return NaN
	endif
	PathInfo folder_path
	Print "•FEDMS_LoadJVFolder(\""+S_path+"\",\""+measurement_persons+"\",\""+measurement_comments+"\",\""+file_format+"\","+num2str(device_area)+","+num2str(mismatch_factor)+")"
	FEDMS_LoadJVFolder(S_path, measurement_persons, measurement_comments, file_format,device_area,mismatch_factor)
End

Function FEDMS_LoadToFDataFolder(isIntensityTest)
	Variable isIntensityTest
	String original_folder = GetDataFolder(1)
	Variable voltage_amplification
	String measurement_persons
	String measurement_comments
	Prompt voltage_amplification, "Enter voltage amplification:"
	Prompt measurement_persons, "Enter the measurement person(s):"
	Prompt measurement_comments, "Enter any additional measurement comments:"
	DoPrompt "Enter Measurement Info", voltage_amplification, measurement_persons, measurement_comments
	if(V_flag==1)
		return NaN
	endif
	NewPath/O/Q folder_path
	if(V_flag==1)
		return NaN
	endif
	String file_list = IndexedFile(folder_path,-1,".dat")
	// Filter list
	Variable i
	for(i=0;i<ItemsInList(file_list,";");i+=1)
		String file_name = StringFromList(i,file_list,";")
		if(StringMatch(file_name,"!electrons*") && StringMatch(file_name,"!holes*"))
			file_list = RemoveListItem(i,file_list)
			i -= 1
		endif
	endfor
	file_list = SortList(file_list,";",16)
	Variable file_count = ItemsInList(file_list)
	String sample_name = FEDMS_ChooseSample()
	if(StringMatch(sample_name,""))
		return NaN
	endif
	SetDataFolder root:FlexEDMS:$(sample_name)
	NewDataFolder/S/O $("Time of Flight")
	// Load data from each file
	Print "Loading "+num2str(ItemsInList(file_list,";"))+" data files..."
	for(i=0;i<ItemsInList(file_list,";");i+=1)
		file_name = StringFromList(i,file_list,";")
		PathInfo folder_path
		String fullpathname = S_path+file_name
		FEDMS_LoadTOFDataFile(isIntensityTest,sample_name,fullpathname,voltage_amplification,measurement_persons,measurement_comments)
	endfor
	SetDataFolder original_folder
End

Function FEDMS_LoadTOFDataFiles()
	String original_folder = GetDataFolder(1)
	Variable voltage_amplification = 0
	String measurement_persons
	String measurement_comments
	Prompt voltage_amplification, "Enter voltage amplification:"
	Prompt measurement_persons, "Enter the measurement person(s):"
	Prompt measurement_comments, "Enter any additional measurement comments:"
	do
		DoPrompt "Enter Measurement Info", voltage_amplification, measurement_persons, measurement_comments
		if(V_flag==1)
			return NaN
		endif
		if(!voltage_amplification>0)
			DoAlert 0, "The voltage amplification value must be greater than 0."
		endif
	while(!voltage_amplification>0)
	String file_list = FEDMS_DoOpenMultiFileDialog()
	if(StringMatch(file_list,""))
		return NaN
	endif
	String sample_name = FEDMS_ChooseSample()
	if(StringMatch(sample_name,""))
		return NaN
	endif
	SetDataFolder root:FlexEDMS:$(sample_name)
	NewDataFolder/S/O $("Time of Flight")
	// Load data from each file
	Print "Loading "+num2str(ItemsInList(file_list,"\r"))+" data files..."
	Variable i
	for(i=0;i<ItemsInList(file_list,"\r");i+=1)
		String fullpathname = StringFromList(i,file_list,"\r")
		FEDMS_LoadTOFDataFile(0,sample_name,fullpathname,voltage_amplification,measurement_persons,measurement_comments)
	endfor
	SetDataFolder original_folder
End

Function FEDMS_LoadTOFDataFile(isIntensityTest,sample_name,fullpathname,amplification,persons,comments)
	Variable isIntensityTest
	String sample_name
	String fullpathname
	Variable amplification
	String persons
	String comments
	String carrier_type
	String measurement_date
	Variable temperature
	Variable bias
	Variable voltage_applied
	Variable intensity
	Variable optical_density
	Variable time_start
	Variable time_step
	String original_folder = GetDataFolder(1)
	SetDataFolder root:FlexEDMS:$(sample_name):$("Time of Flight"):
	// Load file into text wave
	String filename = StringFromList(ItemsInList(fullpathname,":")-1,fullpathname,":")
	LoadWave/A/J/Q/K=2/V={""," $",0,0} fullpathname
	Wave/T wave0
	// Parse filename
	carrier_type = StringFromList(0,filename,"_")
	NewDataFolder/S/O $carrier_type
	if(isIntensityTest)
		NewDataFolder/S/O intensity_test
	endif
	// Parse file header
	Variable j
	Variable data_start_index
	for(j=0;j<numpnts(wave0);j+=1)
		if(!StringMatch(wave0[j],"#*"))
			data_start_index = j
			break
		else
			String header_list = StringFromList(2,wave0[j],"#")
			Variable k
			for(k=0;k<ItemsInList(header_list);k+=1)
				String var_name = StringFromList(0,StringFromList(k,header_list),"=")
				String var_value = StringFromList(1,StringFromList(k,header_list),"=")
				if(StringMatch(var_name,"date"))
					measurement_date = var_value
				endif
				if(StringMatch(var_name,"T"))
					temperature = str2num(var_value)
				endif
				if(StringMatch(var_name,"tstart"))
					time_start = str2num(var_value)
				endif
				if(StringMatch(var_name,"dt"))
					time_step = str2num(var_value)
				endif
				if(StringMatch(var_name,"Vp"))
					voltage_applied = str2num(var_value)
				endif
				if(StringMatch(var_name,"Intensity"))
					intensity = str2num(var_value)
				endif
				if(StringMatch(var_name,"ODvK"))
					optical_density = str2num(var_value)
				endif
				if(StringMatch(var_name,"remark"))
					if(StringMatch(comments,""))
						comments += var_value
					else
						comments += ","+var_value
					endif
				endif
			endfor
		endif
	endfor
	KillWaves wave0
	// Create data folders
	NewDataFolder/S/O $num2istr(temperature)+"K"
	NewDataFolder/S/O $num2str(amplification*voltage_applied)+"V"
	if(isIntensityTest)
		NewDataFolder/S/O $num2str(intensity/10^(optical_density))+"%"
	endif
	// Create global measurement variables and strings
	Variable/G Voltage_amplification = amplification
	String/G Measurement_persons = persons
	String/G Measurement_comments = comments
	// Load raw voltage measurement wave
	LoadWave/A/J/D/Q/K=0/L={0,data_start_index,0,0,0} fullpathname
	Duplicate/O $("wave0") voltage_raw
	SetScale/P x time_start,time_step,"s", voltage_raw
	// Delete the points before time=0
	Variable V_baseline = mean(voltage_raw,0,50)
	voltage_raw -= V_baseline
	DeletePoints 0,x2pnt(voltage_raw,0), voltage_raw
	SetScale/P x 0,time_step,"s", voltage_raw
	// Filter excess data points to obtain log spaced data
	Variable pnts_per_dec = 50
	Variable log_range = log(pnt2x(voltage_raw,numpnts(voltage_raw)-1))-log(pnt2x(voltage_raw,1))
	Variable log_size = round(log_range*pnts_per_dec)
	Make/N=(log_size)/D/O time_log, V_norm
	for(j=0;j<numpnts(time_log);j+=1)
		time_log[j] = 10^(log(pnt2x(voltage_raw,1))+(1/pnts_per_dec)*j)
		V_norm[j] = abs(voltage_raw(time_log[j]))
	endfor
	// Normalize data to peak value
	Variable max_val =  WaveMax(V_norm)
	V_norm /= max_val
	// Cleanup
	KillWaves $("wave0") voltage_raw
	SetDataFolder original_folder
End

Function FEDMS_PlotCapacitanceVsBias(device_name) : Graph
	String device_name
	String frequency_name
	String original_data_folder = GetDataFolder(1)
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[Strlen(device_name)-1]
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance
	Wave/T illuminations
	SetDataFolder :$("CV_"+illuminations[0]):
	Display/N=$("Capacitance_Data_"+device_name) $("C_mu") vs $("voltage_applied")
	Variable i
	for(i=1;i<numpnts(illuminations);i+=1)
		SetDataFolder ::$("CV_"+illuminations[i]):
		AppendToGraph $("C_mu") vs $("voltage_applied")
	endfor
	Execute "FEDMS_GraphStyle()"
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance
	Wave C_sats
	if(WaveExists(C_sats))
		AppendToGraph $("C_sats") vs $("V_sats")
		ModifyGraph mode(C_sats)=3, marker(C_sats)=18, rgb(C_sats)=(0,0,0)
	endif
	TextBox/C/N=text0/F=0/A=LT/X=8/Y=8 device_name
	Label left "Chemical Capacitance (F)"
	Label bottom "Voltage (V)"
	SetAxis/A=2 left
	SetAxis bottom -1,*
	ModifyGraph log(bottom)=0
	SetDataFolder original_data_folder
End

Function FEDMS_PlotCarrierDensityVsBias(device_name) : Graph
	String device_name
	String original_data_folder = GetDataFolder(1)
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[Strlen(device_name)-1]
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance
	Wave/T illuminations
	Variable i
	for(i=0;i<numpnts(illuminations);i+=1)
		SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance:$("CV_"+illuminations[i]):
		if(i==0)
			Display/W=(506.4,240.2,744,408.8)/N=$("Bias_Carrier_Density_"+device_name) $("carrier_density") vs $("voltage_applied")
		else
			AppendToGraph $("carrier_density") vs $("voltage_applied")
		endif
	endfor
	Execute "FEDMS_GraphStyle()"
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance
	AppendToGraph $("n_scs") vs $("V_scs")
	AppendToGraph $("n_mps") vs $("V_mps")
	AppendToGraph $("n_ocs") vs $("V_ocs")
	ModifyGraph mode(n_scs)=3, marker(n_scs)=18, rgb(n_scs)=(0,0,0)
	ModifyGraph mode(n_mps)=3, marker(n_mps)=18, rgb(n_mps)=(0,0,0)
	ModifyGraph mode(n_ocs)=3, marker(n_ocs)=18, rgb(n_ocs)=(0,0,0)
	ModifyGraph log(bottom)=0
	ModifyGraph axOffset(left)=-2.25
	ModifyGraph gridRGB=(13056,13056,13056)
	ModifyGraph gridHair=0
	ModifyGraph ZisZ=1
	TextBox/C/N=text0/F=0/A=MT/X=0/Y=6 device_name
	Label left "Charge Carrier Density (cm\\S-3\\M)"
	Label bottom "Voltage (V)"
	SetAxis left 1e14,*
	SetAxis bottom -1,*
	// Add legend
	i = numpnts(illuminations)-1
	Legend/C/N=text1/J/F=0/A=LT/X=4/Y=4 "\\s(carrier_density#"+num2str(i)+") "+illuminations[i]
	for(i-=1;i>=0;i-=1)
		if(i==0)
			AppendText "\\s(carrier_density) "+illuminations[i];DelayUpdate
		else
			AppendText "\\s(carrier_density#"+num2str(i)+") "+illuminations[i];DelayUpdate
		endif
	endfor
	SetDataFolder original_data_folder
End

Function FEDMS_PlotFFVsTheta(device_list) : Graph
	String device_list
	String original_data_folder = GetDataFolder(1)
	Variable i
	for(i=0;i<ItemsInList(device_list);i+=1)
		String device_name = StringFromList(i,device_list)	
		String substrate_num = RemoveEnding(device_name)
		String device_num = device_name[Strlen(device_name)-1]
		SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):IPDA_Results
		if(i==0)
			Display $("FFs") vs $("theta")
		else
			AppendToGraph $("FFs") vs $("theta")
		endif
	endfor
	Execute "FEDMS_GraphStyle()"
	ShowKBColorizePanel()
	KBColorizeTraces#KBColorTablePopMenuProc("name",0,"YellowHot")
	KillWindow $("KBColorizePanel")
	SetAxis left 0.3,0.75
	ModifyGraph log(left)=0
	SetAxis bottom 0.04,40
	Label left "Fill Factor"
	Label bottom "θ"
	SetDataFolder original_data_folder
End
End

Function FEDMS_PlotKbrVsBias(device_name) : Graph
	String device_name
	String original_data_folder = GetDataFolder(1)
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[Strlen(device_name)-1]
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance
	Wave/T illuminations
	Variable i
	for(i=0;i<numpnts(illuminations);i+=1)
		SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance:$("CV_"+illuminations[i]):
		if(i==0)
			Display/W=(506.4,240.2,744,408.8)/N=$("Bias_Kbr_"+device_name) $("k_br") vs $("voltage_applied")
		else
			AppendToGraph $("k_br") vs $("voltage_applied")
		endif
	endfor
	Execute "FEDMS_GraphStyle()"
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance
	AppendToGraph $("k_mps") vs $("V_mps")
	AppendToGraph $("k_ocs") vs $("V_ocs")
	ModifyGraph mode(k_mps)=3, marker(k_mps)=18, rgb(k_mps)=(0,0,0)
	ModifyGraph mode(k_ocs)=3, marker(k_ocs)=18, rgb(k_ocs)=(0,0,0)
	ModifyGraph log(bottom)=0
	ModifyGraph axOffset(left)=-2.25
	ModifyGraph gridRGB=(13056,13056,13056)
	ModifyGraph gridHair=0
	ModifyGraph ZisZ=1
	TextBox/C/N=text0/F=0/A=MT/X=0/Y=6 device_name
	Label left "Recombination Coefficient, k\Bbr\M (cm\\S3\\Ms\\S-1\\M)"
	Label bottom "Voltage (V)"
	SetAxis/A=2 left
	SetAxis bottom 0,*
	// Add legend
	for(i=0;i<numpnts(illuminations);i+=1)
		if(i==0)
			Legend/C/N=text1/J/F=0/A=LT "\\s(k_br) "+illuminations[i]
		else
			AppendText "\\s(k_br#"+num2str(i)+") "+illuminations[i];DelayUpdate
		endif
	endfor
	SetDataFolder original_data_folder
End

Function FEDMS_PlotMobilityVsBias(device_name) : Graph
	String device_name
	String original_data_folder = GetDataFolder(1)
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[Strlen(device_name)-1]
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance
	Wave/T illuminations
	Variable i
	for(i=0;i<numpnts(illuminations);i+=1)
		SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance:$(illuminations[i]):
		if(i==0)
			Display/W=(506.4,240.2,744,408.8)/N=$("Bias_Mobility_"+device_name) $("mu_eff") vs $("voltage_applied")
		else
			AppendToGraph $("mu_eff") vs $("voltage_applied")
		endif
	endfor
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):Impedance
	Execute "FEDMS_GraphStyle()"
	AppendToGraph $("mu_mps") vs $("V_mps")
	ModifyGraph mode(mu_mps)=3, marker(mu_mps)=18, rgb(mu_mps)=(0,0,0)
	ModifyGraph log(bottom)=0
	ModifyGraph axOffset(left)=-2.25
	ModifyGraph gridRGB=(13056,13056,13056)
	ModifyGraph gridHair=0
	ModifyGraph ZisZ=1
	TextBox/C/N=text0/F=0/A=MT/X=0/Y=6 device_name
	Label left "Charge Carrier Mobility (cm\\S2\\MV\\S-1\\Ms\\S-1\\M)"
	Label bottom "Voltage (V)"
	SetAxis/A=2 left
	SetAxis bottom 0,*
	// Add legend
	for(i=0;i<numpnts(illuminations);i+=1)
		if(i==0)
			Legend/C/N=text1/J/F=0/A=LT "\\s(mu_eff) "+illuminations[i]
		else
			AppendText "\\s(mu_eff#"+num2str(i)+") "+illuminations[i];DelayUpdate
		endif
	endfor
	SetDataFolder original_data_folder
End

Function FEDMS_PlotMobilityVsField(device_name) : Graph
	String device_name
	String original_data_folder = GetDataFolder(1)
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[Strlen(device_name)-1]
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):IPDA_Results
	Wave n_field
	Wave mu_vs_field
	Wave field_sqrt
	Variable i
	for(i=0;i<numpnts(n_field);i+=1)
		if(i==0)
			Display mu_vs_field[][0] vs field_sqrt[][0]
		else
			AppendToGraph mu_vs_field[][i] vs field_sqrt[][i]
		endif
	endfor
	Execute "FEDMS_GraphStyle()"
	SetAxis left 1e-005,0.001
	SetAxis bottom 0,350
	ModifyGraph log(bottom)=0
	ModifyGraph mode=4,marker=19
	TextBox/C/N=text0/F=0/Z=1/A=MT device_name
	Legend/C/N=text1/J/F=0/A=LT "Charge carrier density (cm\S-3\M)"
	for(i=0;i<numpnts(n_field);i+=1)
		if(i==0)
			AppendText "\\s(mu_vs_field)"+num2str(n_field[i])
		else
			AppendText "\\s(mu_vs_field#"+num2str(i)+")"+num2str(n_field[i])
		endif
	endfor
	Label left "Charge Carrier Mobility (cm\\S2\\MV\\S-1\\Ms\\S-1\\M)"
	Label bottom "Electric Field, F\S1/2\M (V/cm)\S1/2\M"
	SetDataFolder original_data_folder
End

Function FEDMS_PlotTOFData(data_path,sample_name,carrier_type,temperature_name,bias_name)
	String data_path
	String sample_name
	String carrier_type
	String temperature_name
	String bias_name
	String original_folder = GetDataFolder(1)
	SetDataFolder data_path
	Wave V_norm
	Wave time_log
	Display V_norm vs time_log
	Execute "FEDMS_GraphStyle()"
	SetAxis left 0.01,2
	WaveStats/Q V_norm
	FindLevel/EDGE=2/R=[V_maxloc]/Q V_norm, 0.01
	SetAxis bottom 1e-7,2*time_log[ceil(V_LevelX)]
	Label left "Normalized Current (a.u.)"
	Label bottom "Time (s)"
	TextBox/C/N=text0/F=0/Z=1/A=MT sample_name+", "+carrier_type+" @ T="+temperature_name+", V="+bias_name
	SetDataFolder original_folder
End

Function FEDMS_PlotTOF_IntensityTest()
	String original_folder = GetDataFolder(1)
	String sample_name = FEDMS_ChooseTOFSample()
	if(StringMatch(sample_name,""))
		return NaN
	endif
	String carrier_type = FEDMS_ChooseTOFCarrierType(sample_name)
	if(StringMatch(carrier_type,""))
		return NaN
	endif
	SetDataFolder root:FlexEDMS:$(sample_name):$("Time of Flight"):$(carrier_type)
	if(!DataFolderExists("intensity_test"))
		DoAlert, 0, "There is no intensity test data for the selected sample and carrier type."
		return NaN
	endif
	SetDataFolder :$("intensity_test")
	String temperature_name = FEDMS_ChooseSubfolder(GetDataFolder(1))
	SetDataFolder :$(temperature_name)
	Wave mobility_geo
	Wave dispersion
	Wave intensity
	Wave voltage_V
	String bias_name = FEDMS_ChooseSubfolder(GetDataFolder(1))
	SetDataFolder :$(bias_name)
	// Build intensity list
	DFREF dfr1 = GetDataFolderDFR()
	Variable N_intensities = CountObjectsDFR(dfr1,4)
	String intensity_list = ""
	Variable i
	for(i=0;i<N_intensities;i+=1)
		intensity_list = AddListItem(GetIndexedObjNameDFR(dfr1,4,i),intensity_list)
	endfor
	intensity_list = SortList(intensity_list,";",16)
	// Plot ToF transients at each intensity
	Variable plot_count = 0
	for(i=0;i<N_intensities;i+=1)
		String intensity_name = StringFromList(i,intensity_list)
		SetDataFolder :$(intensity_name)
		Wave V_norm
		Wave time_log
		if(plot_count==0)
			FEDMS_PlotTOFData(GetDataFolder(1),sample_name,carrier_type,temperature_name,bias_name)
		else
			AppendToGraph V_norm vs time_log
		endif
		SetDataFolder ::
		plot_count += 1
	endfor
	Execute "FEDMS_GraphStyle()"
	KBColorizeTracesLinearHue(0.5, 1, 0)
	// Plot mobility vs intensity and dispersion vs intensity
	SetDataFolder ::
	Variable voltage = str2num(RemoveEnding(bias_name,"V"))
	Variable index_start = -1
	Variable index_end = -1
	for(i=0;i<numpnts(voltage_V);i+=1)
		if(index_start<0 && voltage_V[i]==voltage)
			index_start = i
		endif
		if(!(index_start<0) && voltage_V[i]!=voltage_V[index_start])
			index_end = i-1
			break
		endif
		if(!(index_start<0) && i==numpnts(voltage_V)-1)
			index_end = i
			break
		endif
	endfor
	Display mobility_geo[index_start,index_end] vs intensity[index_start,index_end]
	Execute "FEDMS_GraphStyle()"
	Label left "Charge Carrier Mobility (cm\\S2\\MV\\S-1\\Ms\\S-1\\M)"
	Label bottom "Intensity"
	Display dispersion[index_start,index_end] vs intensity[index_start,index_end]
	Execute "FEDMS_GraphStyle()"
	Label left "Dispersion"
	Label bottom "Intensity"
	// Plot dispersion vs intensity
	SetDataFolder original_folder
End

Function FEDMS_PlotJV(device_name,measurement_type,measurement_name) : Graph
	String device_name
	String measurement_type
	String measurement_name
	String original_data_folder = GetDataFolder(1)
	// Parse the device name
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[strlen(device_name)-1]
	// Check that the device folder exists
	if(!DataFolderExists("root:FlexEDMS:"+substrate_num+":"+device_num))
		Print "Error! Device name not found. Check that the sample and data have been correctly loaded."
		return NaN
	endif
	// Analyze the specified J-V measurement
	FEDMS_AnalyzeJV(device_name,measurement_type,measurement_name)
	// Begin plotting the results
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):$(measurement_type):$(measurement_name):
	Wave J_avg
	if(StringMatch(measurement_name,"*dark*"))
		Display /W=(8.25,46.25,414,339.5) $"J_avg_abs"
	else
		Display  /W=(426.75,46.25,832.5,339.5) J_avg
	endif
	Execute "FEDMS_GraphStyle()"
	ModifyGraph log=0, zero=1, width=360, margin(left)=36, margin(bottom)=32, margin(right)=10, margin(top)=10
	// Custom formatting for dark J-V plots
	if(StringMatch(measurement_name,"*dark*"))
		NVAR R_s, R_sh
		ModifyGraph log(left)=1
		TextBox/C/N=text1/F=0/A=MC ("R_s = "+num2str(R_s)+" Ohm \rR_sh = "+num2str(R_sh)+" Ohm")
	// Custom formatting for light J-V plots
	else
		NVAR J_sc, V_oc, P_max, FF
		SetAxis left *,5
		FindLevel/Q J_avg,5
		SetAxis bottom *,(V_LevelX+0.1)
		TextBox/C/N=text1/F=0/A=MC ("J_sc = "+num2str(J_sc)+" mA/cm^2 \rV_oc = "+num2str(V_oc)+" V \rP_max = "+num2str(P_max)+" mW/cm^2 \rFF = "+num2str(FF))
	endif
	Label left "Current Density (mA cm\\S-2\\M)"
	Label bottom "Voltage (V)"
	TextBox/C/N=text0/F=0/A=MT/X=0.00/Y=6.00 ("Device #"+device_name+" - "+measurement_name)
	SetDataFolder original_data_folder
End

Function FEDMS_PlotJVGUI() : Graph
	String device_name = FEDMS_ChooseDevice("JV")
	if(StringMatch(device_name,""))
		return NaN
	endif
	String measurement_name = FEDMS_ChooseMeasurement(device_name,"JV")
	Print "•FEDMS_PlotJV(\""+device_name+"\",\"JV\",\""+measurement_name+"\")"
	FEDMS_PlotJV(device_name,"JV",measurement_name)
End

Function FEDMS_PlotPhotocurrentData(device_name,[show_fits])
	String device_name
	Variable show_fits
	String original_data_folder = GetDataFolder(1)
	String substrate_num = RemoveEnding(device_name)
	String device_num = device_name[Strlen(device_name)-1]
	// Check that the device folder exists
	if(!DataFolderExists("root:FlexEDMS:"+substrate_num+":"+device_num))
		Print "Error! Device name not found. Check that the sample and data have been correctly loaded."
		return NaN
	endif
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV
	Wave/T illuminations
	Wave J_sat_onsets
	Wave V_eff_onsets
	// Plot first photocurrent curve
	SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV:$(illuminations[0]):
	// Skip plotting if photocurrent data is missing
	if(Exists("photocurrent")!=1)
		Print "Error! No photocurrent data found for that sample. Make sure you run the FEDMS_ first"
		SetDataFolder original_data_folder
		return NaN
	endif
	Display /W=(580,40,1080,365) $("photocurrent") vs $("voltage_effective")
	// Append the rest of the photocurrent curves
	Variable i
	for(i=1;i<numpnts(illuminations);i+=1)
		SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV:$(illuminations[i]):
		AppendToGraph $("photocurrent") vs $("voltage_effective")
	endfor
	Execute "FEDMS_GraphStyle()"
	SetAxis bottom 0.01,*
	SetAxis/A=2 left
	Label left "Photocurrent Density (mA cm\\S-2\\M)"
	Label bottom "Effective Voltage (V)"
	ModifyGraph log=1, mode=0, lsize=2
	TextBox/C/N=text0/F=0/A=MB/X=4/Y=4 (device_name)
	if(show_fits==1)
		for(i=0;i<numpnts(illuminations);i+=1)
			SetDataFolder root:FlexEDMS:$(substrate_num):$(device_num):JV:$(illuminations[i]):
			if(Exists("fit_photocurrrent")==1)
				AppendToGraph $("fit_photocurrent") vs $("voltage_effective")
				ModifyGraph rgb[numpnts(illuminations)+i] = (0,0,0)
			endif
		endfor
		if(WaveExists(J_sat_onsets))
			AppendToGraph J_sat_onsets vs V_eff_onsets
			ModifyGraph mode(J_sat_onsets)=3, marker(J_sat_onsets)=18, rgb(J_sat_onsets)=(0,0,0)
		endif
	endif
	// Add Legend
	i = numpnts(illuminations)-1
	Legend/C/N=text1/J/F=0/A=RB/X=4/Y=4 "\\s(photocurrent#"+num2str(i)+") "+illuminations[i]
	for(i-=1;i>=0;i-=1)
		if(i==0)
			AppendText "\\s(photocurrent) "+illuminations[i];DelayUpdate
		else
			AppendText "\\s(photocurrent#"+num2str(i)+") "+illuminations[i];DelayUpdate
		endif
	endfor
	// Restore working directory
	SetDataFolder original_data_folder
End

Function FEDMS_PlotPhotocurrentDataGUI()
	String device_name = FEDMS_ChooseDevice("JV")
	if(StringMatch(device_name,""))
		return NaN
	endif
	Print "•FEDMS_AnalyzeJVIntensity(\""+device_name+"\")"
	FEDMS_AnalyzeJVIntensity(device_name)
	Print "•FEDMS_PlotPhotocurrentData(\""+device_name+"\",show_fits=1)"
	FEDMS_PlotPhotocurrentData(device_name,show_fits=1)
End

Function FEDMS_PlotTOF_FieldDependences()
	String sample_name = FEDMS_ChooseTOFSample()
	if(StringMatch(sample_name,""))
		return NaN
	endif
	String carrier_type = FEDMS_ChooseTOFCarrierType(sample_name)
	if(StringMatch(carrier_type,""))
		return NaN
	endif
	String original_folder = GetDataFolder(1)
	SetDataFolder root:FlexEDMS:$(sample_name):$("Time of Flight"):$(carrier_type)
	// Build the temperature list
	String temperature_option
	DFREF dfr1 = GetDataFolderDFR()
	Variable N_temps = CountObjectsDFR(dfr1,4)
	String temperature_list = "all"
	Variable i
	for(i=0;i<N_temps;i+=1)
		String temperature_name = GetIndexedObjNameDFR(dfr1,4,i)
		if(StringMatch(temperature_name,"!intensity_test"))
			temperature_list = AddListItem(temperature_name,temperature_list)
		endif
	endfor
	temperature_list = SortList(temperature_list,";",16)
	// Prompt user to choose temperature
	Prompt temperature_option, "Choose a temperature option:" popup, temperature_list
	DoPrompt "Make Selection", temperature_option
	// Loop though all temperatures and plot mobility vs electric_field
	// Invert temperature list to plot highest temp first
	temperature_list = SortList(temperature_list,";",17)
	temperature_list = RemoveEnding(temperature_list,";all")
	Variable plot_count = 0
	for(i=0;i<N_temps;i+=1)
		temperature_name = StringFromList(i,temperature_list)
		// Check for selected temperature option
		if(StringMatch(temperature_option,"all") || StringMatch(temperature_name,temperature_option))
			SetDataFolder :$(temperature_name)
			Wave mobility_geo
			Wave electric_field
			if(plot_count==0)
				Display mobility_geo vs electric_field
			else
				AppendToGraph mobility_geo vs electric_field
			endif
			SetDataFolder ::
			plot_count += 1
		endif
	endfor
	Execute "FEDMS_GraphStyle()"
	ModifyGraph mode=4,marker=19,log(bottom)=0,gaps=0
	Label left "Charge Carrier Mobility (cm\\S2\\MV\\S-1\\Ms\\S-1\\M)"
	Label bottom "Electric Field (V/m)"
	TextBox/C/N=text0/F=0/A=LB sample_name+"\r"+carrier_type
	// Loop though all temperatures and plot dispersion vs electric_field
	plot_count = 0
	for(i=0;i<N_temps;i+=1)
		temperature_name = StringFromList(i,temperature_list)
		// Check for selected temperature option
		if(StringMatch(temperature_option,"all") || StringMatch(temperature_name,temperature_option))
			SetDataFolder :$(temperature_name)
			Wave dispersion
			if(plot_count==0)
				Display dispersion vs electric_field
			else
				AppendToGraph dispersion vs electric_field
			endif
			SetDataFolder ::
			plot_count += 1
		endif
	endfor
	Execute "FEDMS_GraphStyle()"
	KBColorizeTracesLinearHue(0.5, 1, 0)
	ModifyGraph mode=4,marker=19,log=0,gaps=0
	Label left "Dispersion"
	Label bottom "Electric Field (V/m)"
	TextBox/C/N=text0/F=0/A=RB sample_name+"\r"+carrier_type
	SetDataFolder original_folder
End

Function FEDMS_PlotTOF_TempDependence()
	String sample_name = FEDMS_ChooseTOFSample()
	if(StringMatch(sample_name,""))
		return NaN
	endif
	String carrier_type = FEDMS_ChooseTOFCarrierType(sample_name)
	if(StringMatch(carrier_type,""))
		return NaN
	endif
	FEDMS_AnalyzeTOF_TempDependence(sample_name,carrier_type)
	String original_folder = GetDataFolder(1)
	SetDataFolder root:FlexEDMS:$(sample_name):$("Time of Flight"):$(carrier_type)
	// Plot mobility vs T
	Wave mobility_avg
	Wave mobility_stdev
	Wave temperature_K
	Wave dispersion_avg
	Wave dispersion_stdev
	Display mobility_avg vs temperature_K
	Execute "FEDMS_GraphStyle()"
	ModifyGraph log(bottom)=0
	ModifyGraph mode=4,marker=19,rgb=(0,0,65535)
	ErrorBars mobility_avg Y,wave=(mobility_stdev,mobility_stdev)
	Label left "Charge Carrier Mobility (cm\\S2\\MV\\S-1\\Ms\\S-1\\M)"
	Label bottom "Temperature (K)"
	TextBox/C/N=text0/F=0/A=LB sample_name+"\r"+carrier_type
	// Plot dispersion vs T
	Display dispersion_avg vs temperature_K
	Execute "FEDMS_GraphStyle()"
	ModifyGraph log=0
	ModifyGraph mode=4,marker=19,rgb=(0,0,65535)
	ErrorBars dispersion_avg Y,wave=(dispersion_stdev,dispersion_stdev)
	Label left "Dispersion"
	Label bottom "Temperature (K)"
	TextBox/C/N=text0/F=0/A=LB sample_name+"\r"+carrier_type
End

Macro FEDMS_GraphStyle() : GraphStyle
	// Axes format
	ModifyGraph log=1,tick=2,mirror=1,logLabel=2,standoff=0,logTicks=3,logHTrip(left)=100,logLTrip(left)=0.01
	// Figure size
	ModifyGraph width=400,height={Aspect,0.7}
	// Figure margins
	ModifyGraph margin(left)=41,margin(bottom)=36,margin(right)=10,margin(top)=10
	// Line and marker style
	ModifyGraph lsize=2
	// Line Colors
	ShowKBColorizePanel()
	KBColorizeTraces#KBColorTablePopMenuProc("name",0,"Rainbow")
	KillWindow $("KBColorizePanel")
End
