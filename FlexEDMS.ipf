#pragma rtGlobals=1	// Use modern global access method.
#pragma IgorVersion = 6.3 // Minimum Igor version required
#pragma version = 0.1.3-alpha 

// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the FlexEDMS project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The FlexEDMS project can be found on Github at https://github.com/MikeHeiber/FlexEDMS

#include <KBColorizeTraces>

Menu "FlexEDMS"
	"Create New Sample", /Q, FEDMS_CreateNewSample()
	"Edit Sample Info", /Q, FEDMS_EditSampleInfo()
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
	ModifyGraph expand=1.5
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
	Make/N=(N_temps)/D/O temperature_K, mobility_avg, mobility_stdev
	String temperature_list = ""
	Variable i
	for(i=0;i<N_temps;i+=1)
		temperature_list = AddListItem(GetIndexedObjNameDFR(dfr1,4,i),temperature_list)
	endfor
	temperature_list = SortList(temperature_list,";",16)
	Make/D/N=5/O mobility_temp
	for(i=0;i<N_temps;i+=1)
		String temperature_name = StringFromList(i,temperature_list)
		SetDataFolder :$(temperature_name)
		Duplicate/O dispersion dispersion_temp
		Wave mobility_geo
		Variable j
		for(j=0;j<numpnts(mobility_temp);j+=1)
			Variable k
			Variable current_val = 1
			Variable index = 0
			for(k=0;k<numpnts(dispersion_temp);k+=1)
				if(dispersion_temp[k]>0 && dispersion_temp[k]<current_val)
					current_val = dispersion_temp[k]
					index = k
				endif
			endfor
			mobility_temp[j] = mobility_geo[index]
			dispersion_temp[index] = 1
		endfor
		temperature_K[i] = str2num(RemoveEnding(temperature_name,"K"))
		WaveStats/Q mobility_temp
		mobility_avg[i] = V_avg
		mobility_stdev[i] = V_sdev
		KillWaves dispersion_temp
		SetDataFolder ::
	endfor
	KillWaves mobility_temp
	SetDataFolder original_folder
End

Function FEDMS_ButtonProc_CancelTOF(ba) : ButtonControl
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

Function FEDMS_ButtonProc_DoneTOF(ba) : ButtonControl
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

Function FEDMS_ButtonProc_FitTOF(ba) : ButtonControl
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
			KillWaves difference W_coef W_sigma log_trace_wave log_times
			// Enter results into table
			String measurement_folder = GetDataFolder(1)
			Variable isIntensityTest = StringMatch(StringFromList(5,measurement_folder,":"),"intensity_test")
			NVAR Measurement_index
			if(isIntensityTest)
				SetDataFolder ::::::
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

Function FEDMS_ButtonProc_RejectTOF(ba) : ButtonControl
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

Function/S FEDMS_CreateNewSample()
	String original_folder = GetDataFolder(1)
	NewDataFolder/O/S root:FlexEDMS
	// Build the current sample list
	String sample_list = ""
	DFREF dfr1 = GetDataFolderDFR()
	Variable N_samples = CountObjectsDFR(dfr1,4)
	Variable i
	for(i=0;i<N_samples;i+=1)
		sample_list = AddListItem(GetIndexedObjNameDFR(dfr1,4,i),sample_list)
	endfor
	String sample_name = ""
	String fab_date
	String fab_person
	String sample_comp
	Variable sample_thickness = 0
	String comments
	Prompt sample_name, "Enter the sample name:"
	Prompt fab_date, "Enter the fabrication date:"
	Prompt fab_person, "Enter the fabrication person:"
	Prompt sample_comp, "Enter the composition of the semiconductor layer:"
	Prompt sample_thickness, "Enter the semiconductor layer thickness (m):"
	Prompt comments, "Enter any additional comments and notes:"
	Variable input_needed = 1
	do
		DoPrompt "Enter Sample Info",sample_name, fab_date, fab_person, sample_comp, sample_thickness, comments
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
		return ""
	endif
	if(!DataFolderExists(sample_name))
		NewDataFolder $(sample_name)
	endif
	SetDataFolder :$(sample_name)
	String/G Fabrication_date = fab_date
	String/G Fabrication_person = fab_person
	Variable/G Thickness_m = sample_thickness
	String/G Sample_composition = sample_comp
	String/G Sample_comments = comments
	SetDataFolder original_folder
	return sample_name
End

Function FEDMS_CreateNewSample2()
	String sample_type_list = FEDMS_LoadSampleTypes()
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

Function FEDMS_EditSampleInfo()
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
		if(StringMatch(folder_name,"Sample*"))
			sample_list = AddListItem(folder_name,sample_list)
		endif
	endfor
	// Prompt user to choose the sample
	String sample_name
	Prompt sample_name, "Choose the sample name:", popup, sample_list
	DoPrompt "Make Selections",sample_name
	// Create sample edit form
	SetDataFolder :$(sample_name)
	SVAR Fabrication_date
	SVAR Fabrication_person
	NVAR Thickness_m
	SVAR Sample_composition
	SVAR Sample_comments
	String fab_date = Fabrication_date
	String fab_person = Fabrication_person
	String sample_comp = Sample_composition
	Variable sample_thickness = Thickness_m
	String comments = Sample_comments
	Prompt sample_name, "Enter the sample name:"
	Prompt fab_date, "Enter the fabrication date:"
	Prompt fab_person, "Enter the fabrication person:"
	Prompt sample_comp, "Enter the composition of the semiconductor layer:"
	Prompt sample_thickness, "Enter the semiconductor layer thickness (m):"
	Prompt comments, "Enter any additional comments and notes:"
	DoPrompt "Enter Sample Info",sample_name, fab_date, fab_person, sample_comp, sample_thickness, comments
	// Update global variables and strings
	Fabrication_date = fab_date
	Fabrication_person = fab_person
	Thickness_m = sample_thickness
	Sample_composition = sample_comp
	Sample_comments = comments
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

Function/S FEDMS_LoadSampleTypes()
	
End

Function FEDMS_LoadToFDataFolder(isIntensityTest)
	Variable isIntensityTest
	String original_folder = GetDataFolder(1)
	Variable voltage_amplification
	String measurement_person
	String measurement_comments
	Prompt voltage_amplification, "Enter voltage amplification:"
	Prompt measurement_person, "Enter the measurement person:"
	Prompt measurement_comments, "Enter any additional measurement comments:"
	DoPrompt "Enter Measurement Info", voltage_amplification, measurement_person, measurement_comments
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
		FEDMS_LoadTOFDataFile(isIntensityTest,sample_name,fullpathname,voltage_amplification,measurement_person,measurement_comments)
	endfor
	SetDataFolder original_folder
End

Function FEDMS_LoadTOFDataFiles()
	String original_folder = GetDataFolder(1)
	Variable voltage_amplification = 0
	String measurement_person
	String measurement_comments
	Prompt voltage_amplification, "Enter voltage amplification:"
	Prompt measurement_person, "Enter the measurement person:"
	Prompt measurement_comments, "Enter any additional measurement comments:"
	do
		DoPrompt "Enter Measurement Info", voltage_amplification, measurement_person, measurement_comments
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
		FEDMS_LoadTOFDataFile(0,sample_name,fullpathname,voltage_amplification,measurement_person,measurement_comments)
	endfor
	SetDataFolder original_folder
End

Function FEDMS_LoadTOFDataFile(isIntensityTest,sample_name,fullpathname,amplification,person,comments)
	Variable isIntensityTest
	String sample_name
	String fullpathname
	Variable amplification
	String person
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
	String/G Measurement_person = person
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
	Variable plot_count = 0
	for(i=0;i<N_intensities;i+=1)
		String intensity_name = StringFromList(i,intensity_list)
		SetDataFolder :$(intensity_name)
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
	SetDataFolder original_folder
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
		temperature_list = AddListItem(GetIndexedObjNameDFR(dfr1,4,i),temperature_list)
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
		String temperature_name = StringFromList(i,temperature_list)
		// Check for selected temperature option
		if(StringMatch(temperature_option,"all") || StringMatch(temperature_name,temperature_option))
			SetDataFolder :$(temperature_name)
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
	KBColorizeTracesLinearHue(0.5, 1, 0)
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
	Display mobility_avg vs temperature_K
	Execute "FEDMS_GraphStyle()"
	ModifyGraph mode=4,marker=19,rgb=(0,0,65535)
	ErrorBars mobility_avg Y,wave=(mobility_stdev,mobility_stdev)
	Label left "Charge Carrier Mobility (cm\\S2\\MV\\S-1\\Ms\\S-1\\M)"
	Label bottom "Temperature (K)"
	TextBox/C/N=text0/F=0/A=LB sample_name+"\r"+carrier_type
End

Macro FEDMS_GraphStyle() : GraphStyle
	// Axes format
	ModifyGraph log=1,tick=2,mirror=1,logLabel=2,standoff=0,logTicks=3,logHTrip(left)=100,logLTrip(left)=0.01
	// Figure size
	ModifyGraph width=415,height={Aspect,0.7}
	// Figure margins
	ModifyGraph margin(left)=43,margin(bottom)=36
End