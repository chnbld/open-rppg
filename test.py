import rppg
import numpy as np
model = rppg.Model()

# User-provided anthropometric data for BMI calculation
# NOTE: BMI cannot be calculated from video/physiological signals - requires actual measurements
USER_HEIGHT_CM = None  # Set your height in centimeters (e.g., 175)
USER_WEIGHT_KG = None  # Set your weight in kilograms (e.g., 70)

video_path = "../test2.mp4"  # Set path to your MP4 file
result = model.process_video(video_path)
if result:
    print(f"Heart Rate: {result['hr']} BPM")
    
    # Calculate Cardiac Stress Index from HRV metrics
    if result.get('hrv'):
        hrv = result['hrv']
        stress_components = {}
        
        # LF/HF ratio - higher values indicate higher stress (sympathetic dominance)
        # Normal range: 0.5-2.0, stress threshold: >2.0
        lf_hf_ratio = hrv.get('LF/HF', None)
        if lf_hf_ratio is not None:
            # Normalize LF/HF (0-1 scale, where 0.5 = 0 stress, 2.0+ = high stress)
            lf_hf_stress = min(1.0, max(0.0, (lf_hf_ratio - 0.5) / 1.5))
            stress_components['lf_hf'] = lf_hf_stress
        
        # RMSSD - lower values indicate higher stress
        # Normal range: 20-60ms, stress threshold: <20ms
        rmssd = hrv.get('rmssd', None)
        if rmssd is not None:
            # Normalize RMSSD (inverted: lower RMSSD = higher stress)
            rmssd_stress = min(1.0, max(0.0, 1.0 - (rmssd / 60.0)))
            stress_components['rmssd'] = rmssd_stress
        
        # SDNN - lower values indicate higher stress
        # Normal range: 20-60ms, stress threshold: <20ms
        sdnn = hrv.get('sdnn', None)
        if sdnn is not None:
            # Normalize SDNN (inverted: lower SDNN = higher stress)
            sdnn_stress = min(1.0, max(0.0, 1.0 - (sdnn / 60.0)))
            stress_components['sdnn'] = sdnn_stress
        
        # Calculate composite stress index (0-100 scale)
        if stress_components:
            # Weighted average: LF/HF gets more weight (0.5), RMSSD and SDNN each 0.25
            weights = {'lf_hf': 0.5, 'rmssd': 0.25, 'sdnn': 0.25}
            # Calculate weighted sum
            weighted_sum = sum(stress_components.get(k, 0) * weights.get(k, 0) 
                             for k in stress_components.keys())
            # Calculate total weight for available components
            total_weight = sum(weights.get(k, 0) for k in stress_components.keys())
            stress_index = weighted_sum / total_weight if total_weight > 0 else 0
            cardiac_stress = stress_index * 100  # Convert to 0-100 scale
            
            # Classify stress level
            if cardiac_stress < 30:
                stress_level = "Low"
            elif cardiac_stress < 50:
                stress_level = "Moderate"
            elif cardiac_stress < 70:
                stress_level = "High"
            else:
                stress_level = "Very High"
            
            print(f"Cardiac Stress Index: {cardiac_stress:.1f}/100 ({stress_level})")
            result['cardiac_stress'] = float(cardiac_stress)
            result['cardiac_stress_level'] = stress_level
            result['cardiac_stress_components'] = stress_components
    
    # Display Heart Rate Variability (HRV) Metrics
    if result.get('hrv'):
        hrv_data = result['hrv']
        print(f"\n📊 Heart Rate Variability (HRV) Metrics:")
        print(f"{'─' * 60}")
        
        # Time Domain Metrics
        if 'bpm' in hrv_data:
            print(f"  Heart Rate (Peak Method): {hrv_data['bpm']:.2f} BPM")
        if 'ibi' in hrv_data:
            print(f"  Inter-Beat Interval (IBI): {hrv_data['ibi']:.2f} ms")
        if 'sdnn' in hrv_data:
            print(f"  SDNN (Std Dev of NN intervals): {hrv_data['sdnn']:.2f} ms")
        if 'sdsd' in hrv_data:
            print(f"  SDSD (Std Dev of Successive Differences): {hrv_data['sdsd']:.2f} ms")
        if 'rmssd' in hrv_data:
            print(f"  RMSSD (Root Mean Square of Successive Differences): {hrv_data['rmssd']:.2f} ms")
        if 'pnn20' in hrv_data:
            print(f"  pNN20 (Proportion of NN > 20ms): {hrv_data['pnn20']:.4f}")
        if 'pnn50' in hrv_data:
            print(f"  pNN50 (Proportion of NN > 50ms): {hrv_data['pnn50']:.4f}")
        if 'hr_mad' in hrv_data:
            print(f"  HR MAD (Heart Rate Median Absolute Deviation): {hrv_data['hr_mad']:.2f} BPM")
        
        # Poincaré Plot Metrics
        if 'sd1' in hrv_data:
            print(f"  SD1 (Short-term variability): {hrv_data['sd1']:.2f} ms")
        if 'sd2' in hrv_data:
            print(f"  SD2 (Long-term variability): {hrv_data['sd2']:.2f} ms")
        if 's' in hrv_data:
            print(f"  S (Poincaré Plot Area): {hrv_data['s']:.2f}")
        if 'sd1/sd2' in hrv_data:
            print(f"  SD1/SD2 Ratio: {hrv_data['sd1/sd2']:.4f}")
        
        # Frequency Domain Metrics
        if 'VLF' in hrv_data:
            print(f"  VLF (Very Low Frequency Power): {hrv_data['VLF']:.4f}")
        if 'LF' in hrv_data:
            print(f"  LF (Low Frequency Power): {hrv_data['LF']:.4f}")
        if 'HF' in hrv_data:
            print(f"  HF (High Frequency Power): {hrv_data['HF']:.4f}")
        if 'LF/HF' in hrv_data:
            print(f"  LF/HF Ratio: {hrv_data['LF/HF']:.4f}")
        if 'TP' in hrv_data:
            print(f"  TP (Total Power): {hrv_data['TP']:.4f}")
        if 'breathingrate' in hrv_data:
            breathing_rate = hrv_data['breathingrate']
            if breathing_rate:
                br_bpm = breathing_rate * 60 if breathing_rate < 1 else breathing_rate
                print(f"  Breathing Rate: {br_bpm:.2f} breaths/min")
        
        print(f"{'─' * 60}")
    
    # Calculate Cardiac Workload
    if result.get('hr') is not None:
        hr = result['hr']
        
        # Method 1: Heart Rate Reserve (HRR) - percentage of max heart rate capacity used
        # Estimated Max HR = 220 - age (using average age of 40 if unknown)
        # For simplicity, we'll use a standard max HR estimate
        estimated_max_hr = 190  # Conservative estimate for general population
        resting_hr = 60  # Average resting heart rate
        
        # Calculate Heart Rate Reserve percentage
        if estimated_max_hr > resting_hr:
            hrr_percentage = ((hr - resting_hr) / (estimated_max_hr - resting_hr)) * 100
            hrr_percentage = max(0, min(100, hrr_percentage))  # Clamp to 0-100
        else:
            hrr_percentage = 0
        
        # Method 2: Cardiac Workload Index based on HR zones and HRV
        # Define heart rate zones (relative to max HR)
        workload_factors = {
            'resting': (0, 50),      # 0-50% max HR
            'light': (50, 60),       # 50-60% max HR
            'moderate': (60, 70),    # 60-70% max HR
            'vigorous': (70, 85),    # 70-85% max HR
            'maximum': (85, 100)     # 85-100% max HR
        }
        
        hr_percent_of_max = (hr / estimated_max_hr) * 100
        
        # Determine HR zone
        if hr_percent_of_max < 50:
            hr_zone = 'resting'
            zone_workload = 20
        elif hr_percent_of_max < 60:
            hr_zone = 'light'
            zone_workload = 40
        elif hr_percent_of_max < 70:
            hr_zone = 'moderate'
            zone_workload = 60
        elif hr_percent_of_max < 85:
            hr_zone = 'vigorous'
            zone_workload = 80
        else:
            hr_zone = 'maximum'
            zone_workload = 100
        
        # Adjust workload based on HRV if available (lower HRV = higher workload)
        hrv_adjustment = 0
        if result.get('hrv'):
            hrv_data = result['hrv']
            # Lower RMSSD or SDNN indicates higher cardiac workload
            if 'rmssd' in hrv_data and hrv_data['rmssd'] is not None:
                # Normalize: RMSSD < 30ms increases workload
                if hrv_data['rmssd'] < 30:
                    hrv_adjustment = 10  # Increase workload
                elif hrv_data['rmssd'] > 50:
                    hrv_adjustment = -5  # Decrease workload
            
            # High LF/HF ratio also increases workload
            if 'LF/HF' in hrv_data and hrv_data['LF/HF'] is not None:
                if hrv_data['LF/HF'] > 2.0:
                    hrv_adjustment += 5
        
        cardiac_workload = zone_workload + hrv_adjustment
        cardiac_workload = max(0, min(100, cardiac_workload))  # Clamp to 0-100
        
        # Classify workload level
        if cardiac_workload < 30:
            workload_level = "Very Light"
        elif cardiac_workload < 50:
            workload_level = "Light"
        elif cardiac_workload < 70:
            workload_level = "Moderate"
        elif cardiac_workload < 85:
            workload_level = "High"
        else:
            workload_level = "Very High"
        
        print(f"Cardiac Workload Index: {cardiac_workload:.1f}/100 ({workload_level})")
        print(f"  Heart Rate Reserve: {hrr_percentage:.1f}%")
        print(f"  HR Zone: {hr_zone.capitalize()} ({hr_percent_of_max:.1f}% of max HR)")
        
        result['cardiac_workload'] = float(cardiac_workload)
        result['cardiac_workload_level'] = workload_level
        result['heart_rate_reserve'] = float(hrr_percentage)
        result['hr_zone'] = hr_zone
        result['hr_percent_of_max'] = float(hr_percent_of_max)
    
    # Calculate Body Mass Index (BMI) if height and weight are provided
    if USER_HEIGHT_CM is not None and USER_WEIGHT_KG is not None:
        height_m = USER_HEIGHT_CM / 100.0  # Convert cm to meters
        bmi = USER_WEIGHT_KG / (height_m ** 2)
        
        # Classify BMI category
        if bmi < 18.5:
            bmi_category = "Underweight"
        elif bmi < 25.0:
            bmi_category = "Normal"
        elif bmi < 30.0:
            bmi_category = "Overweight"
        else:
            bmi_category = "Obese"
        
        print(f"Body Mass Index (BMI): {bmi:.1f} ({bmi_category})")
        print(f"  Height: {USER_HEIGHT_CM} cm")
        print(f"  Weight: {USER_WEIGHT_KG} kg")
        
        result['bmi'] = float(bmi)
        result['bmi_category'] = bmi_category
        result['height_cm'] = float(USER_HEIGHT_CM)
        result['weight_kg'] = float(USER_WEIGHT_KG)
    # else:
    #     print("⚠️  BMI not calculated: Please set USER_HEIGHT_CM and USER_WEIGHT_KG in the script")
    
    # Extract and display Respiratory Rate from HRV data
    if result.get('hrv') and 'breathingrate' in result['hrv']:
        breathing_rate = result['hrv']['breathingrate']
        if breathing_rate:
            # Check if value is in Hz (breaths/sec) or bpm (breaths/min)
            # If value < 1, assume it's in Hz; if >= 5, assume it's already in bpm
            if breathing_rate < 1:
                respiratory_rate_bpm = breathing_rate * 60  # Convert Hz to bpm
            else:
                respiratory_rate_bpm = breathing_rate  # Already in bpm
            print(f"Respiratory Rate: {respiratory_rate_bpm:.2f} breaths/min")
            result['respiratory_rate'] = respiratory_rate_bpm
    else:
        # Calculate respiratory rate from BVP signal if not in HRV
        print("Respiratory rate not found in HRV, calculating from BVP signal...")
    
    # Get Blood Volume Pulse (BVP) data
    bvp, timestamps = model.bvp()
    if len(bvp) > 0:
        # Calculate BVP summary statistics
        bvp_mean = float(np.mean(bvp))
        bvp_std = float(np.std(bvp))
        bvp_amplitude = float(np.max(bvp) - np.min(bvp))  # Peak-to-peak amplitude
        bvp_peak = float(np.max(np.abs(bvp)))  # Maximum absolute value
        
        print(f"\nBVP Summary:")
        print(f"  Mean: {bvp_mean:.4f}")
        print(f"  Std Dev: {bvp_std:.4f}")
        print(f"  Amplitude (peak-to-peak): {bvp_amplitude:.4f}")
        print(f"  Peak Value: {bvp_peak:.4f}")
        print(f"  Data Points: {len(bvp)}")
        
        # Add summary values to result (keeping full array for detailed analysis)
        result['bvp_mean'] = bvp_mean
        result['bvp_std'] = bvp_std
        result['bvp_amplitude'] = bvp_amplitude
        result['bvp_peak'] = bvp_peak
        result['bvp'] = bvp.tolist() if hasattr(bvp, 'tolist') else list(bvp)  # Full signal
        result['bvp_timestamps'] = timestamps.tolist() if hasattr(timestamps, 'tolist') else list(timestamps)
        
        # Calculate respiratory rate from BVP if not already calculated
        if 'respiratory_rate' not in result:
            try:
                from scipy.signal import welch
                # Respiratory rate is typically 0.15-0.4 Hz (9-24 breaths/min)
                f, Pxx = welch(bvp, fs=model.fps, nperseg=min(len(bvp), 256), nfft=4096)
                respiratory_range = (f >= 0.15) & (f <= 0.4)
                if np.any(respiratory_range):
                    peak_idx = np.argmax(Pxx[respiratory_range])
                    peak_freq = f[respiratory_range][peak_idx]
                    respiratory_rate_bpm = peak_freq * 60
                    print(f"Respiratory Rate (calculated from BVP): {respiratory_rate_bpm:.2f} breaths/min")
                    result['respiratory_rate'] = float(respiratory_rate_bpm)
            except Exception as e:
                print(f"Could not calculate respiratory rate: {e}")
        
        # Estimate Blood Pressure from BVP and HR/HRV
        # try:
        #     from scipy.signal import find_peaks
            
        #     # Analyze pulse wave characteristics
        #     # Find peaks in BVP signal (systolic peaks)
        #     peaks, properties = find_peaks(bvp, distance=int(model.fps * 0.5), prominence=0.1)
            
        #     if len(peaks) > 1:
        #         # Calculate pulse characteristics
        #         pulse_amplitudes = bvp[peaks]
        #         mean_amplitude = np.mean(pulse_amplitudes)
        #         std_amplitude = np.std(pulse_amplitudes)
                
        #         # Calculate pulse intervals (RR intervals in seconds)
        #         if len(peaks) > 2:
        #             pulse_intervals = np.diff(timestamps[peaks])
        #             mean_rr = np.mean(pulse_intervals)
        #         else:
        #             mean_rr = 60.0 / result.get('hr', 70)  # Estimate from HR
                
        #         # Method 1: Pulse Amplitude Index (PAI) - correlates with BP
        #         # Higher pulse amplitude variation may indicate higher BP
        #         amplitude_variation = std_amplitude / (mean_amplitude + 1e-8) if mean_amplitude > 0 else 0
                
        #         # Method 2: HR-based estimation (limited accuracy)
        #         hr = result.get('hr', 70)
        #         # Base BP estimates (mmHg) - using population averages adjusted by HR
        #         base_sbp = 120  # Average systolic
        #         base_dbp = 80  # Average diastolic
                
        #         # Adjust based on HR (higher HR often correlates with higher BP during activity)
        #         hr_adjustment = (hr - 70) * 0.3  # ~0.3 mmHg per BPM above 70
        #         sbp_estimate = base_sbp + hr_adjustment
        #         dbp_estimate = base_dbp + hr_adjustment * 0.6  # Diastolic changes less
                
        #         # Method 3: HRV-based adjustments
        #         hrv_adjustment_sbp = 0
        #         hrv_adjustment_dbp = 0
        #         if result.get('hrv'):
        #             hrv_data = result['hrv']
        #             # Lower HRV (higher stress) often correlates with higher BP
        #             if 'sdnn' in hrv_data and hrv_data['sdnn'] is not None:
        #                 # SDNN < 30ms suggests higher BP
        #                 if hrv_data['sdnn'] < 30:
        #                     hrv_adjustment_sbp = 10
        #                     hrv_adjustment_dbp = 5
        #                 elif hrv_data['sdnn'] > 60:
        #                     hrv_adjustment_sbp = -5
        #                     hrv_adjustment_dbp = -3
                    
        #             # High LF/HF ratio (stress) may correlate with elevated BP
        #             if 'LF/HF' in hrv_data and hrv_data['LF/HF'] is not None:
        #                 if hrv_data['LF/HF'] > 2.0:
        #                     hrv_adjustment_sbp += 5
        #                     hrv_adjustment_dbp += 3
                
        #         # Method 4: Pulse amplitude variation adjustment
        #         # Higher variation may indicate arterial stiffness (higher BP)
        #         amplitude_adjustment = min(15, max(-10, amplitude_variation * 20))
                
        #         # Final estimates
        #         sbp_estimate = sbp_estimate + hrv_adjustment_sbp + amplitude_adjustment
        #         dbp_estimate = dbp_estimate + hrv_adjustment_dbp + amplitude_adjustment * 0.5
                
        #         # Clamp to reasonable ranges
        #         sbp_estimate = max(90, min(180, sbp_estimate))
        #         dbp_estimate = max(50, min(120, dbp_estimate))
                
        #         # Calculate Mean Arterial Pressure (MAP)
        #         map_estimate = dbp_estimate + (sbp_estimate - dbp_estimate) / 3
                
        #         # Calculate Pulse Pressure
        #         pulse_pressure = sbp_estimate - dbp_estimate
                
        #         # Classify BP category
        #         if sbp_estimate < 120 and dbp_estimate < 80:
        #             bp_category = "Normal"
        #         elif sbp_estimate < 130 and dbp_estimate < 80:
        #             bp_category = "Elevated"
        #         elif sbp_estimate < 140 or dbp_estimate < 90:
        #             bp_category = "High Stage 1"
        #         elif sbp_estimate < 180 or dbp_estimate < 120:
        #             bp_category = "High Stage 2"
        #         else:
        #             bp_category = "Hypertensive Crisis"
                
        #         print(f"\n⚠️  Blood Pressure Estimation (UNCALIBRATED - Approximate Only):")
        #         print(f"  Systolic BP: {sbp_estimate:.0f} mmHg")
        #         print(f"  Diastolic BP: {dbp_estimate:.0f} mmHg")
        #         print(f"  MAP: {map_estimate:.0f} mmHg")
        #         print(f"  Pulse Pressure: {pulse_pressure:.0f} mmHg")
        #         print(f"  Category: {bp_category}")
        #         print(f"  ⚠️  Note: These are estimates. For accurate BP, use a calibrated device.")
                
        #         result['blood_pressure'] = {
        #             'systolic': float(sbp_estimate),
        #             'diastolic': float(dbp_estimate),
        #             'map': float(map_estimate),
        #             'pulse_pressure': float(pulse_pressure),
        #             'category': bp_category,
        #             'method': 'pulse_wave_analysis_hr_hrv_estimated',
        #             'accuracy_warning': True
        #         }
        #     else:
        #         print("⚠️  Insufficient BVP signal quality for blood pressure estimation")
                
        # except Exception as e:
        #     print(f"⚠️  Could not estimate blood pressure: {e}")
        
        #print(f"BVP added to result")

# with model.video_capture(0):          # Connect to your webcam
#     while True:
#         result = model.hr(start=-15)  # Get heart rate from last 15 seconds
#         if result:
#             print(f"Heart Rate: {result['hr']} BPM")
#         time.sleep(1)