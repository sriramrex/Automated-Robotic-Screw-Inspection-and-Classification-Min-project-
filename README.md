# Automated-Robotic-Screw-Inspection-and-Classification-Min-project-

# 🤖 Automated Robotic Screw Inspection and Classification

An **industrial automation project** that integrates **robotics, artificial vision, and data-driven classification** to automatically **inspect, measure, classify, and sort screws** using two industrial robots and machine vision.

---

## 📌 Project Overview

This project implements a **fully automated quality inspection system** for screws, combining:

* **Robotic manipulation** (ABB IRB 140 & KUKA KR3)
* **ROS-based coordination**
* **Industrial cameras & structured lighting**
* **MATLAB-based image processing**
* **Geometric feature extraction & classification**

The system performs:

* Screw head inspection
* Screw thread measurement
* Defect detection
* Screw classification & sorting

---

## 🏗️ System Architecture

### 🔧 Hardware Components

* **ABB IRB 140** robot

  * Picks screws from feeder
  * Places them on screw holder
  * Controls conveyor & metrology station

* **KUKA KR3 R540** robot

  * Equipped with a camera
  * Moves to predefined camera poses
  * Captures screw head images

* **Cameras**

  * MatrixVision BlueCougar **2 MP** (Sorting station)
  * MatrixVision BlueCougar **5 MP** (Metrology station)

* **Optics & Illumination**

  * Telecentric lens (metrology)
  * Front light + backlight illumination

* **Motorized screw holder** for 360° inspection

---

### 🧠 Software Stack

| Layer             | Technology                       |
| ----------------- | -------------------------------- |
| Robot Control     | ABB RAPID, ROS (MoveIt)          |
| Communication     | ROS topics, Socket communication |
| Vision Processing | MATLAB                           |
| Motion Planning   | Forward & Inverse Kinematics     |
| Classification    | Feature-based decision logic     |

---

## 🔁 System Workflow

```text
Feeder → ABB Robot → Screw Holder → KUKA Camera Inspection
        ↓                              ↓
     Conveyor → Metrology Station → Measurement & Classification
        ↓
   Sorting & Palletizing
```

---

## 🦾 Robot Modeling & Motion Planning

### ✔ Forward Kinematics

* Computes end-effector pose from joint angles
* Homogeneous transformation matrices used
* Required for accurate camera positioning

### ✔ Inverse Kinematics

* Analytical solution for first 3 joints
* Rotation-matrix-based solution for wrist joints
* Newton–Raphson method for numerical refinement

### ✔ Body Jacobian

* Maps joint velocities to end-effector velocities
* Used for:

  * Velocity control
  * Manipulability analysis

### ✔ Trajectory Planning

* **Trapezoidal velocity profiles**
* Smooth acceleration & deceleration
* Point-to-point and via-point trajectories

---

## 🎥 Vision-Based Screw Inspection

### 🔍 Measured Parameters

* Screw head diameter / width
* Total screw length
* Thread angle
* Thread pitch
* Pitch diameter

---

## 🧪 Metrology Station – Thread Inspection (MATLAB)

### Processing Pipeline

1. Image acquisition
2. Background subtraction
3. Otsu thresholding (binarization)
4. Morphological operations
5. Noise removal (median filtering)
6. Hough Circle detection
7. Hough Line transform
8. Thread angle extraction
9. Pitch & depth calculation
10. Scaling (pixel → mm)

📏 **Scale factor:** ~0.03 mm/pixel

---

## 🧠 Sorting Station – Screw Head Inspection

### Camera Calibration

* MATLAB Camera Calibration App
* Intrinsic & extrinsic parameter estimation
* Homography-based rectification

### Head Measurement

* Hough Circle Transform → round heads
* Hough Line Transform → hex heads
* Dimension conversion to millimeters

---

## 📡 ROS Communication Architecture

### ROS Nodes

* Main controller node (Python)
* KUKA motion control node
* Image acquisition nodes
* MATLAB processing interface
* Classification & analytics node

### Communication Methods

* ROS Publisher / Subscriber
* Socket communication (ABB ↔ ROS)

---

## 📊 Screw Classification

Based on extracted features:

* Metric vs non-metric thread
* Damaged vs non-damaged
* Screw type identification

Results are sent to ABB robot for **sorting & palletizing**.

---

## ⚠️ Known Limitations

| Issue                | Cause                    | Proposed Solution            |
| -------------------- | ------------------------ | ---------------------------- |
| Poor binarization    | Similar background color | Different fixture color      |
| Limited height range | Telecentric lens         | Zoom lens with calibration   |
| Shiny / dark heads   | Lighting issues          | Diffuse + darkfield lighting |

---

## 🚀 Future Improvements

* Single robot with **camera + gripper** on same end-effector
* Motorized focus lens
* Machine learning-based classification
* Faster inspection cycle via parallel processing

---

## 📁 Project Structure (Suggested)

```text
├── robotics/
│   ├── kuka_ros/
│   ├── abb_rapid/
├── vision/
│   ├── metrology_station/
│   ├── sorting_station/
├── matlab/
│   ├── image_processing/
│   ├── calibration/
├── docs/
│   ├── diagrams/
│   ├── report.pdf
├── README.md
```

---

## 🧾 References

* Industrial Robotics – ABB & KUKA
* MATLAB Image Processing Toolbox
* ROS Industrial & MoveIt
* ISO Metric Thread Standards

---

## 👨‍💻 Author

**Automated Robotic Screw Inspection and Classification**
Academic & Industrial Robotics Project

---

⭐ *If you use this project, consider citing or referencing the work.*
