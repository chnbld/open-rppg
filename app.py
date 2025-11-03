import os
import json
import numpy as np
from flask import Flask, render_template, request, jsonify, send_from_directory
from werkzeug.utils import secure_filename
import rppg

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'uploads'
app.config['MAX_CONTENT_LENGTH'] = 500 * 1024 * 1024  # 500MB max file size
app.config['ALLOWED_EXTENSIONS'] = {'mp4', 'avi', 'mov', 'mkv', 'webm'}

# Create upload folder if it doesn't exist
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in app.config['ALLOWED_EXTENSIONS']

def calculate_cardiac_stress(result):
    """Calculate cardiac stress index from HRV metrics"""
    if not result.get('hrv'):
        return result
    
    hrv = result['hrv']
    stress_components = {}
    
    # LF/HF ratio - higher values indicate higher stress
    lf_hf_ratio = hrv.get('LF/HF', None)
    if lf_hf_ratio is not None:
        lf_hf_stress = min(1.0, max(0.0, (lf_hf_ratio - 0.5) / 1.5))
        stress_components['lf_hf'] = lf_hf_stress
    
    # RMSSD - lower values indicate higher stress
    rmssd = hrv.get('rmssd', None)
    if rmssd is not None:
        rmssd_stress = min(1.0, max(0.0, 1.0 - (rmssd / 60.0)))
        stress_components['rmssd'] = rmssd_stress
    
    # SDNN - lower values indicate higher stress
    sdnn = hrv.get('sdnn', None)
    if sdnn is not None:
        sdnn_stress = min(1.0, max(0.0, 1.0 - (sdnn / 60.0)))
        stress_components['sdnn'] = sdnn_stress
    
    # Calculate composite stress index (0-100 scale)
    if stress_components:
        weights = {'lf_hf': 0.5, 'rmssd': 0.25, 'sdnn': 0.25}
        weighted_sum = sum(stress_components.get(k, 0) * weights.get(k, 0) 
                         for k in stress_components.keys())
        total_weight = sum(weights.get(k, 0) for k in stress_components.keys())
        stress_index = weighted_sum / total_weight if total_weight > 0 else 0
        cardiac_stress = stress_index * 100
        
        # Classify stress level
        if cardiac_stress < 30:
            stress_level = "Low"
        elif cardiac_stress < 50:
            stress_level = "Moderate"
        elif cardiac_stress < 70:
            stress_level = "High"
        else:
            stress_level = "Very High"
        
        result['cardiac_stress'] = float(cardiac_stress)
        result['cardiac_stress_level'] = stress_level
        result['cardiac_stress_components'] = stress_components
    
    return result

def calculate_cardiac_workload(result):
    """Calculate cardiac workload index"""
    if result.get('hr') is None:
        return result
    
    hr = result['hr']
    estimated_max_hr = 190
    resting_hr = 60
    
    # Calculate Heart Rate Reserve percentage
    if estimated_max_hr > resting_hr:
        hrr_percentage = ((hr - resting_hr) / (estimated_max_hr - resting_hr)) * 100
        hrr_percentage = max(0, min(100, hrr_percentage))
    else:
        hrr_percentage = 0
    
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
    
    # Adjust workload based on HRV
    hrv_adjustment = 0
    if result.get('hrv'):
        hrv_data = result['hrv']
        if 'rmssd' in hrv_data and hrv_data['rmssd'] is not None:
            if hrv_data['rmssd'] < 30:
                hrv_adjustment = 10
            elif hrv_data['rmssd'] > 50:
                hrv_adjustment = -5
        
        if 'LF/HF' in hrv_data and hrv_data['LF/HF'] is not None:
            if hrv_data['LF/HF'] > 2.0:
                hrv_adjustment += 5
    
    cardiac_workload = zone_workload + hrv_adjustment
    cardiac_workload = max(0, min(100, cardiac_workload))
    
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
    
    result['cardiac_workload'] = float(cardiac_workload)
    result['cardiac_workload_level'] = workload_level
    result['heart_rate_reserve'] = float(hrr_percentage)
    result['hr_zone'] = hr_zone
    result['hr_percent_of_max'] = float(hr_percent_of_max)
    
    return result

def extract_respiratory_rate(result):
    """Extract respiratory rate from HRV data (most accurate)"""
    if result.get('hrv') and 'breathingrate' in result['hrv']:
        breathing_rate = result['hrv']['breathingrate']
        if breathing_rate is not None:
            # Check if it's in Hz (value < 1) or already in bpm
            if breathing_rate < 1:
                respiratory_rate_bpm = breathing_rate * 60
            else:
                respiratory_rate_bpm = breathing_rate
            
            # Validate range (normal breathing: 8-25 breaths/min)
            if 5 <= respiratory_rate_bpm <= 30:
                result['respiratory_rate'] = float(respiratory_rate_bpm)
                result['respiratory_rate_source'] = 'hrv'
    
    return result

def calculate_respiratory_rate_from_bvp(bvp, fps, result):
    """Calculate respiratory rate from BVP signal (fallback method)"""
    try:
        from scipy.signal import welch, find_peaks
        
        # Use longer window for better frequency resolution
        nperseg = min(len(bvp), 512)
        nfft = max(4096, nperseg * 4)
        
        f, Pxx = welch(bvp, fs=fps, nperseg=nperseg, nfft=nfft)
        
        # Rest mode: normal breathing is 0.15-0.35 Hz (9-21 breaths/min)
        # More restrictive range for rest mode to avoid noise
        respiratory_range = (f >= 0.15) & (f <= 0.35)
        
        if np.any(respiratory_range):
            # Get power in respiratory range
            respiratory_power = Pxx[respiratory_range]
            respiratory_freqs = f[respiratory_range]
            
            # Find peaks in the respiratory range (more robust than just max)
            peak_indices, peak_properties = find_peaks(
                respiratory_power,
                prominence=np.max(respiratory_power) * 0.1,  # At least 10% of max power
                distance=5  # Minimum distance between peaks
            )
            
            if len(peak_indices) > 0:
                # Use the peak with highest prominence
                prominences = peak_properties['prominences']
                best_peak_idx = peak_indices[np.argmax(prominences)]
                peak_freq = respiratory_freqs[best_peak_idx]
            else:
                # Fallback to maximum power
                peak_idx = np.argmax(respiratory_power)
                peak_freq = respiratory_freqs[peak_idx]
            
            respiratory_rate_bpm = peak_freq * 60
            
            # Validate range
            if 8 <= respiratory_rate_bpm <= 25:
                result['respiratory_rate'] = float(respiratory_rate_bpm)
                result['respiratory_rate_source'] = 'bvp_spectral'
            elif 5 <= respiratory_rate_bpm <= 30:
                # Allow slightly wider range but mark as potentially less accurate
                result['respiratory_rate'] = float(respiratory_rate_bpm)
                result['respiratory_rate_source'] = 'bvp_spectral'
                result['respiratory_rate_warning'] = 'Rate outside normal rest range'
    except Exception as e:
        pass
    
    return result

def add_bvp_summary(result, model):
    """Add BVP summary statistics and full signal to result"""
    try:
        bvp, timestamps = model.bvp()
        if len(bvp) > 0:
            bvp_mean = float(np.mean(bvp))
            bvp_std = float(np.std(bvp))
            bvp_amplitude = float(np.max(bvp) - np.min(bvp))
            bvp_peak = float(np.max(np.abs(bvp)))
            
            result['bvp_mean'] = bvp_mean
            result['bvp_std'] = bvp_std
            result['bvp_amplitude'] = bvp_amplitude
            result['bvp_peak'] = bvp_peak
            
            # Include full BVP signal and timestamps for visualization
            result['bvp'] = bvp.tolist() if hasattr(bvp, 'tolist') else list(bvp)
            result['bvp_timestamps'] = timestamps.tolist() if hasattr(timestamps, 'tolist') else list(timestamps)
            
            # Calculate respiratory rate if not already extracted from HRV
            if 'respiratory_rate' not in result:
                result = calculate_respiratory_rate_from_bvp(bvp, model.fps, result)
    except Exception as e:
        pass
    
    return result

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload_file():
    if 'video' not in request.files:
        return jsonify({'error': 'No video file provided'}), 400
    
    file = request.files['video']
    if file.filename == '':
        return jsonify({'error': 'No file selected'}), 400
    
    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        file.save(filepath)
        return jsonify({'success': True, 'filename': filename})
    
    return jsonify({'error': 'Invalid file type'}), 400

@app.route('/process', methods=['POST'])
def process_video():
    try:
        data = request.json
        filename = data.get('filename')
        
        if not filename:
            return jsonify({'error': 'No filename provided'}), 400
        
        filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        
        if not os.path.exists(filepath):
            return jsonify({'error': 'File not found'}), 404
        
        # Initialize model and process video
        model = rppg.Model()
        result = model.process_video(filepath)
        
        if not result:
            return jsonify({'error': 'Failed to process video'}), 500
        
        # Add additional calculations
        result = calculate_cardiac_stress(result)
        result = calculate_cardiac_workload(result)
        
        # Extract respiratory rate from HRV first (most accurate), then BVP as fallback
        result = extract_respiratory_rate(result)
        result = add_bvp_summary(result, model)
        
        # Convert numpy arrays to lists for JSON serialization
        def convert_to_serializable(obj):
            if isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, dict):
                return {key: convert_to_serializable(value) for key, value in obj.items()}
            elif isinstance(obj, list):
                return [convert_to_serializable(item) for item in obj]
            return obj
        
        result = convert_to_serializable(result)
        
        return jsonify({'success': True, 'result': result})
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename)

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)

