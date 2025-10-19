# 🫀 open-rppg
An easy-to-use rPPG inference toolbox

## ⚙️ Installation 
Python >= 3.9 & <= 3.12
```bash
pip install open-rppg
```
## 🧪 Import and Use 
```python
import rppg
model  = rppg.Model()
result = model.process_video("your_video.mp4")
```

## 📊 Result Example 
```python
{'hr': 100.45004500450045,        # Heart Rate (FFT Method)
 'SQI': 0.8665749931341046,       # Signal Quality
 'hrv':{
    'bpm': 103.6194862200475,     # Heart Rate (Peak Method)
    'ibi': 579.0416666666666,     # Inter-Beat Interval
    'sdnn': 54.76628055757589,    # Standard Deviation of NN intervals
    'sdsd': 30.674133962201175,   # Standard Deviation of Successive Differences
    'rmssd': 46.25344260031846,   # Root Mean Square of Successive Differences
    'pnn20': 0.5714285714285714,  # Proportion of NN50 > 20ms
    'pnn50': 0.2857142857142857,  # Proportion of NN50 > 50ms
    'hr_mad': 8.333333333333314,  # Heart Rate Median Absolute Deviation
    'sd1': 29.276576197229755,    # Short-term variability
    'sd2': 59.75143144804733,     # Long-term variability
    's': 5495.642490576809,       # Poincaré Plot Area
    'sd1/sd2': 0.489972800445545, # SD1/SD2 Ratio
    'breathingrate': 0.21607605877268798,
    'VLF': 0.09521664913596516,   # Very Low Frequency Power
    'TP': 2.056694418632364,      # Total Power
    'HF': 1.2267116642737315,     # High Frequency Power
    'LF': 0.7347661052226675,     # Low Frequency Power
    'LF/HF': 0.5989721355243509   # LF/HF Ratio
  },
 'latency': 0.0}                  # Real-Time Latency
```

## 🕒 Real-Time Mode 
```python
import time
model = rppg.Model()

with model.video_capture(0):          # Connect to your webcam
    while True:
        result = model.hr(start=-15)  # Get heart rate from last 15 seconds
        if result:
            print(f"Heart Rate: {result['hr']} BPM")
        time.sleep(1)
```

## 🖼️ Real-Time Frame Preview

```python
for frame, box in model.preview:     # Current RGB frame and detection box
    x, y  = box                      
    face  = frame[x[0]:x[1], y[0]:y[1]]
```

## 💓 Get BVP Wave 
```python
bvp, ts        = model.bvp()         # BVP with timestampes
raw_bvp, ts    = model.bvp(raw=True) # Unfiltered BVP
```

## ⌛ Time Slice 
```python
now       = model.now                      # Video duration or current time
bvp, ts   = model.bvp(start=10, end=20)    # BVP slice from 10 to 20 seconds
bvp, ts   = model.bvp(start=-15)           # The last 15-second slice
hr        = model.hr(start=-15)            # HR of the last 15 seconds 
```

## 🧠 Model Selection 
```python
print(rppg.supported_models) # ['ME-chunk.rlap', 'ME-flow.rlap', .......]
model = rppg.Model('RhythmMamba.rlap') # RhythmMamba trained on rlap
```
## 🧰 Pretrained Models 
| Model | Training Set | Description | Paper |
|-|-|-|-| 
|ME-chunk|PURE RLAP|rPPG based on state-space model|[2025](https://doi.org/10.48550/arXiv.2504.01774)|
|ME-flow|PURE RLAP|ME in low-latency real-time mode|[2025](https://doi.org/10.48550/arXiv.2504.01774)| 
|PhysMamba|PURE RLAP|Mamba with fast-slow network|[2024](https://doi.org/10.48550/arXiv.2409.12031)|
|RhythmMamba|PURE RLAP|Mamba with 1D FFT|[2025](https://doi.org/10.1609/aaai.v39i10.33204)|
|PhysFormer|PURE RLAP|Transformer with central diff conv|[2022](https://openaccess.thecvf.com/content/CVPR2022/papers/Yu_PhysFormer_Facial_Video-Based_Physiological_Measurement_With_Temporal_Difference_Transformer_CVPR_2022_paper.pdf)| 
|TSCAN|PURE RLAP|Conv attention with temporal shift|[2020](https://papers.nips.cc/paper/2020/file/e1228be46de6a0234ac22ded31417bc7-Paper.pdf)|
|EfficientPhys|PURE RLAP|TSCAN with self attention|[2022](https://openaccess.thecvf.com/content/WACV2023/papers/Liu_EfficientPhys_Enabling_Simple_Fast_and_Accurate_Camera-Based_Cardiac_Measurement_WACV_2023_paper.pdf)|
|PhysNet|PURE RLAP|3D CNN encoder-decoder network|[2019](https://bmvc2019.org/wp-content/uploads/papers/0186-paper.pdf)| 

## ⚡ Use CUDA 
Install JAX with CUDA (Linux only).
```bash
pip install jax[cuda]
```

## 📜 Licensing Notice

This repository includes source code and tools released under the [MIT License](LICENSE).

However, **pretrained models and model configurations** provided in this repository are the intellectual property of their respective authors and are licensed under the terms specified by the **original papers**. Please refer to the respective publications and their repositories for license details before using these models in your work.

We do **not** claim ownership or rights to redistribute third-party models unless explicitly stated.

## 📚 Citation

```bibtex
@article{yu2019remote,
  title={Remote photoplethysmograph signal measurement from facial videos using spatio-temporal networks},
  author={Yu, Zitong and Li, Xiaobai and Zhao, Guoying},
  journal={arXiv preprint arXiv:1905.02419},
  year={2019}
}

@article{liu2020multi,
  title={Multi-task temporal shift attention networks for on-device contactless vitals measurement},
  author={Liu, Xin and Fromm, Josh and Patel, Shwetak and McDuff, Daniel},
  journal={Advances in Neural Information Processing Systems},
  volume={33},
  pages={19400--19411},
  year={2020}
}

@inproceedings{liu2023efficientphys,
  title={Efficientphys: Enabling simple, fast and accurate camera-based cardiac measurement},
  author={Liu, Xin and Hill, Brian and Jiang, Ziheng and Patel, Shwetak and McDuff, Daniel},
  booktitle={Proceedings of the IEEE/CVF winter conference on applications of computer vision},
  pages={5008--5017},
  year={2023}
}

@inproceedings{yu2022physformer,
  title={Physformer: Facial video-based physiological measurement with temporal difference transformer},
  author={Yu, Zitong and Shen, Yuming and Shi, Jingang and Zhao, Hengshuang and Torr, Philip HS and Zhao, Guoying},
  booktitle={Proceedings of the IEEE/CVF conference on computer vision and pattern recognition},
  pages={4186--4196},
  year={2022}
}

@inproceedings{luo2024physmamba,
  title={PhysMamba: Efficient Remote Physiological Measurement with SlowFast Temporal Difference Mamba},
  author={Luo, Chaoqi and Xie, Yiping and Yu, Zitong},
  booktitle={Chinese Conference on Biometric Recognition},
  pages={248--259},
  year={2024},
  organization={Springer}
}

@inproceedings{zou2025rhythmmamba,
  title={RhythmMamba: Fast, Lightweight, and Accurate Remote Physiological Measurement},
  author={Zou, Bochao and Guo, Zizheng and Hu, Xiaocheng and Ma, Huimin},
  booktitle={Proceedings of the AAAI Conference on Artificial Intelligence},
  volume={39},
  number={10},
  pages={11077--11085},
  year={2025}
}

@article{wang2025memory,
  title={Memory-efficient Low-latency Remote Photoplethysmography through Temporal-Spatial State Space Duality},
  author={Wang, Kegang and Tang, Jiankai and Fan, Yuxuan and Ji, Jiatong and Shi, Yuanchun and Wang, Yuntao},
  journal={arXiv preprint arXiv:2504.01774},
  year={2025}
}
```
