# **ALICE Experiment Data Analysis**
## **Description**

This project is part of a research effort to analyze data from the ALICE detector complex at the Large Hadron Collider (LHC). The project is written in C++ using the ROOT mathematical package, which is widely used in high-energy physics, and includes a large number of tools for statistical processing, visualization and presentation of data. Accordingly, we can distinguish two main parts of the project: the file LHC17d.root contains the initial data, their processing is realized in the file LHC17d.cpp. The result of the project is several important experimental dependencies, which are used in further parts of the study.

## Сompilation

**Attention!** ROOT package must be installed for compilation! Detailed installation instructions on the official website at the link: [***ROOT package installation instructions***](https://root.cern/install/).
In order to compile a project, you need to perform a number of actions:
1. Copy the repository
 
```bash
git clone https://github.com/DoctorVoid3000/My_Projects
 ```
2. Specify the path to the project:

```bash
cd My_Projects/C++/ALICE_Experiment_Data_Analysis/src
```
3. Compile the project using the built-in ROOT compiler:

```bash
root LHC17d.cpp
 ```

