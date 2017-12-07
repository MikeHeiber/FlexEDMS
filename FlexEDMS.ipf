#pragma rtGlobals=1	// Use modern global access method.
#pragma IgorVersion = 6.3 // Minimum Igor version required
#pragma version = 0.1-alpha 

// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the FlexEDMS project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The FlexEDMS project can be found on Github at https://github.com/MikeHeiber/FlexEDMS

Menu "FlexEDMS"
	"Create New Sample", CreateNewSample()
	"Edit Sample Info", EditSampleInfo()
	"Analyze ToF", AnalyzeTOF()
	"Load ToF Data", LoadTOFData()
End

Function AnalyzeTOF()
	String original_folder = GetDataFolder(1)
	String sample_name
	String carrier_type
	// Build the sample list
	SetDataFolder root:FlexEDMS
	String sample_list = ""
	DFREF dfr1 = GetDataFolderDFR()
	Variable N_samples = CountObjectsDFR(dfr1,4)
	String folder_name
	Variable i
	for(i=0;i<N_samples;i+=1)
		folder_name = GetIndexedObjNameDFR(dfr1,4,i)
		if(StringMatch(folder_name,"*"))
			SetDataFolder :$(folder_name)
			if(DataFolderExists("Time of Flight"))
				sample_list = AddListItem(folder_name,sample_list)
			endif
			SetDataFolder ::
		endif
	endfor
	// Prompt user to choose the sample
	Prompt sample_name, "Choose the sample name:", popup, sample_list
	DoPrompt "Make Selections",sample_name
	// User cancelled operation
	if(V_flag==1)
		SetDataFolder original_folder
		return NaN
	endif
	// Build the carrier type list
	SetDataFolder :$(sample_name):$("Time of Flight")
	String carrier_type_list = ""
	DFREF dfr2 = GetDataFolderDFR()
	Variable N_types = CountObjectsDFR(dfr2,4)
	for(i=0;i<N_types;i+=1)
		folder_name = GetIndexedObjNameDFR(dfr2,4,i)
		carrier_type_list = AddListItem(folder_name,carrier_type_list)
	endfor
	Prompt carrier_type, "Choose the carrier type:", popup, carrier_type_list
	DoPrompt "Make Selections",carrier_type
	if(V_flag==1)
		SetDataFolder original_folder
		return NaN
	endif
	SetDataFolder :$(carrier_type)
	// Create needed data waves and variables
	Variable/G V_status = 1
	Wave/Z temperature_K
	if(!WaveExists(temperature_K))
		Make/N=1/D $"temperature_K"
		Wave temperature_K
	endif
	Wave/Z voltage_V
	if(!WaveExists(voltage_V))
		Make/N=1/D $"voltage_V"
		Wave voltage_V
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
	DFREF dfr1 = GetDataFolderDFR()
	Variable N_temps = CountObjectsDFR(dfr1,4)
	String temp_list = ""
	for(i=0;i<N_temps;i+=1)
		temp_list = AddListItem(GetIndexedObjNameDFR(dfr1,4,i),temp_list)
	endfor
	temp_list = SortList(temp_list,";",16)
	Variable measurement_count = 0
	for(i=0;i<N_temps;i+=1)
		String temp_name = StringFromList(i,temp_list)
		SetDataFolder :$(temp_name)
		DFREF dfr2 = GetDataFolderDFR()
		Variable N_biases = CountObjectsDFR(dfr2,4)
		Variable j
		for(j=0;j<N_biases;j+=1)
			String bias_name = GetIndexedObjNameDFR(dfr2,4,j)
			SetDataFolder :$(bias_name)
			Variable/G Measurement_index = measurement_count
			temperature_K[Measurement_index] = {str2num(StringFromList(0,temp_name,"K"))}
			voltage_V[Measurement_index] = {str2num(StringFromList(0,bias_name,"V"))}
			AnalyzeTOFData(sample_name,temp_name,carrier_type,bias_name)
			if(V_status==-1)
				KillVariables V_status
				return NaN
			endif
			measurement_count+=1
			SetDataFolder ::
		endfor
		SetDataFolder ::
	endfor
	KillVariables V_status
	SetDataFolder original_folder
End

Function AnalyzeTOFData(sample_name,temp_name,carrier_type,bias_name)
	String sample_name
	String temp_name
	String carrier_type
	String bias_name
	String original_folder = GetDataFolder(1)
	SetDataFolder root:FlexEDMS:$(sample_name):$("Time of Flight"):$(carrier_type):$(temp_name):$(bias_name)
	PlotTOFData(sample_name,carrier_type,temp_name,bias_name)
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
	DrawText 0,0,"\\Z10Move cursors to pre- and post-transit regions and click Analyze."
	Button button0 title="Analyze",proc=FitTOFButtonProc,size={60,25},fSize=14
	Button button1 title="Done",proc=DoneTOFButtonProc,size={60,25},fSize=14
	Button button2 title="Cancel",proc=CancelTOFButtonProc,size={60,25},fSize=14
	PauseForUser $graph_name
	SetDataFolder original_folder
End

Function CancelTOFButtonProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba
	switch( ba.eventCode )
		case 2: // mouse up
			// click code here
			SetDataFolder :::
			NVAR V_status
			V_status = -1
			KillWindow $WinName(0,1)
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function/S CreateNewSample()
	String original_folder = GetDataFolder(1)
	String sample_name
	String fab_date
	String fab_person
	String sample_comp
	Variable sample_thickness
	String comments
	Prompt sample_name, "Enter the sample name:"
	Prompt fab_date, "Enter the fabrication date:"
	Prompt fab_person, "Enter the fabrication person:"
	Prompt sample_comp, "Enter the composition of the semiconductor layer:"
	Prompt sample_thickness, "Enter the semiconductor layer thickness (m):"
	Prompt comments, "Enter any additional comments and notes:"
	DoPrompt "Enter Sample Info",sample_name, fab_date, fab_person, sample_comp, sample_thickness, comments
	// User cancelled operation
	if(V_flag==1)
		SetDataFolder original_folder
		return ""
	endif
	SetDataFolder root:FlexEDMS
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

Function CreateNewSample2()
	LoadSampleTypes()
End

Function DoneTOFButtonProc(ba) : ButtonControl
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

Function/S DoOpenMultiFileDialog()
	Variable refNum
	String message = "Select one or more files"
	String outputPaths
	String fileFilters = "Data Files (*.txt,*.dat,*.csv):.txt,.dat,.csv;"
	fileFilters += "All Files:.*;"
	Open /D /R /MULT=1 /F=fileFilters /M=message refNum
	outputPaths = S_fileName
	if (strlen(outputPaths) == 0)
		Print "Cancelled"
	endif
	return outputPaths		// Will be empty if user canceled
End

Function EditSampleInfo()
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

Function FilterMobilityData(sample_name,dispersion_limit)
	String sample_name
	Variable dispersion_limit
	String original_folder = GetDataFolder(1)
	SetDataFolder root:FlexEDMS:$(sample_name)
	Wave dispersion
	Wave temperature_K
	Wave mobility_geo
	Make/N=1/D/O mobility_filtered, temp_filtered, dispersion_filtered
	Variable i
	Variable count = 0
	for(i=0;i<numpnts(mobility_geo);i+=1)
		if(dispersion[i]<dispersion_limit)
			dispersion_filtered[count] = {dispersion[i]}
			temp_filtered[count] = {temperature_K[i]}
			mobility_filtered[count] = {mobility_geo[i]}
			count += 1
		endif
	endfor
End

Function FitTOFButtonProc(ba) : ButtonControl
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
			CurveFit/N/Q line log_trace_wave[pcsr(A)-10,pcsr(A)+10] /X=log_times /D
			Wave W_coef
			fit_pre_transit = 10^(W_coef[0]+W_coef[1]*log_times)
			Duplicate/O W_coef W_coef_pre
			if(ItemsInList(trace_list)==1)
				AppendToGraph fit_pre_transit vs times
				ModifyGraph lstyle(fit_pre_transit)=3,rgb(fit_pre_transit)=(0,0,0)
			endif
			CurveFit/N/Q line log_trace_wave[pcsr(B)-10,pcsr(B)+10] /X=log_times /D
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
			Variable dispersion_w
			Duplicate/O fit_pre_transit difference
			difference = fit_pre_transit-fit_post_transit
			FindLevel/Q difference, 0
			transit_time_geo = times(V_LevelX)
			difference = fit_pre_transit_half-fit_post_transit
			FindLevel/Q difference, 0
			transit_time_half = times(V_LevelX)
			dispersion_w = (transit_time_half-transit_time_geo)/transit_time_half
			TextBox/C/N=text1/F=0/Z=1/A=LB "t\\Bgeo\\M = "+num2str(transit_time_geo)+"\rt\\B1/2\\M = "+num2str(transit_time_half)+"\rw = "+num2str(dispersion_w)
			// Clean up from fitting
			KillWaves difference W_coef W_sigma log_trace_wave log_times 
			// Enter results into table
			String measurement_folder = GetDataFolder(1)
			NVAR Measurement_index
			SetDataFolder :::::
			NVAR Thickness_m
			SetDataFolder measurement_folder
			SetDataFolder :::
			Wave voltage_V
			Wave transit_time_geo_s
			Wave transit_time_half_s
			Wave dispersion
			Wave mobility_geo
			transit_time_geo_s[Measurement_index] = {transit_time_geo}
			transit_time_half_s[Measurement_index] = {transit_time_half}
			dispersion[Measurement_index] = {dispersion_w}
			mobility_geo[Measurement_index] = {(Thickness_m*1e2)^2/(abs(voltage_V[Measurement_index])*transit_time_geo)}
			SetDataFolder measurement_folder
			break
		case -1: // control being killed
			break
	endswitch
	return 0
End

Function LoadSampleTypes()
	
End

Function LoadTOFData()
	String original_folder = GetDataFolder(1)
	Variable voltage_amplification
	String measurement_person
	String measurement_comments
	Prompt voltage_amplification, "Enter voltage amplification:"
	Prompt measurement_person, "Enter the measurement person:"
	Prompt measurement_comments, "Enter any additional measurement comments:"
	DoPrompt "Enter Measurement Info", voltage_amplification, measurement_person, measurement_comments
	String file_list = DoOpenMultiFileDialog()
	String sample_name
	String folder_path
	SetDataFolder root:FlexEDMS
	String sample_list = ""
	DFREF dfr1 = GetDataFolderDFR()
	Variable N_samples = CountObjectsDFR(dfr1,4)
	String folder_name
	sample_list = AddListItem("Add New Sample",sample_list)
	Variable i
	for(i=0;i<N_samples;i+=1)
		folder_name = GetIndexedObjNameDFR(dfr1,4,i)
		if(StringMatch(folder_name,"Sample*"))
			sample_list = AddListItem(folder_name,sample_list)
		endif
	endfor
	Prompt sample_name, "Choose the sample:", popup, sample_list
	DoPrompt "Make Selections",sample_name
	if(StringMatch(sample_name,"Add New Sample"))
		sample_name = CreateNewSample()
	endif
	SetDataFolder root:FlexEDMS:$(sample_name)
	NewDataFolder/S/O $("Time of Flight")
	// Load data from each file
	Print "Loading "+num2str(ItemsInList(file_list,"\r"))+" data files..."
	for(i=0;i<ItemsInList(file_list,"\r");i+=1)
		String fullpathname = StringFromList(i,file_list,"\r")
		LoadTOFDataFile(sample_name,fullpathname,voltage_amplification,measurement_person,measurement_comments)
	endfor
	SetDataFolder original_folder
End

Function LoadTOFDataFile(sample_name,fullpathname,amplification,person,comments)
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
	Variable time_start
	Variable time_step
	String original_folder = GetDataFolder(1)
	SetDataFolder root:FlexEDMS:$(sample_name):$("Time of Flight"):
	// Load file into text wave
	String filename = StringFromList(ItemsInList(fullpathname,":")-1,fullpathname,":")
	LoadWave/A/J/Q/K=2 fullpathname
	Wave/T wave0
	// Parse filename
	carrier_type = StringFromList(0,filename,"_")
	NewDataFolder/S/O $carrier_type
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
End

Function PlotTOFData(sample_name,carrier_type,temp_name,bias_name)
	String sample_name
	String carrier_type
	String temp_name
	String bias_name
	String original_folder = GetDataFolder(1)
	SetDataFolder root:FlexEDMS:$(sample_name):$("Time of Flight"):$(carrier_type):$(temp_name):$(bias_name)
	String time_waves = WaveList("t*",";","")
	Wave times = $(StringFromList(0,time_waves))
	String voltage_waves = WaveList("*V_*",";","")
	Wave voltage = $(StringFromList(0,voltage_waves))
	Display voltage vs times;DelayUpdate
	ModifyGraph log=1,tick=2,mirror=1,logLabel=1,standoff=0
	SetAxis left 0.01,2
	SetAxis bottom 1e-7,4e-4
	Label left "Normalized Current (a.u.)"
	Label bottom "Time (s)"
	TextBox/C/N=text0/F=0/Z=1/A=MT carrier_type+" @ T="+temp_name+", V="+bias_name
	ModifyGraph expand=1.5
	SetDataFolder original_folder
End

