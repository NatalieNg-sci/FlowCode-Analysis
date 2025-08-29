# FlowCode-Analysis
For the downstream analysis of FlowCode data after it has been decoded using this - https://github.com/DrCytometer/FlowcodeDecoder.
Based on a FlowCode analysis script written by Dr Oliver Burton. Sections annotated as such were contributed by Dr Vaclav Gergelits. The other sections were written by me.

The main script is FlowCode_Analysis_Natalie, and is used to process a type of T cells you have previously gated, such as CD8s or Tregs. 
Then to compare the T cell types you can use the Compare_T_Cell_Types script. 
And if you had two seperate groups of micein your FlowCode experiment that you want to identify differences in genetic dependancies in, you can use the Compare_2_Groups script. I used it for comparing COPD vs healthy mice, but it can be modified for other purposes easily.
