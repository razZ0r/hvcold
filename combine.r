combined_data <- cbind(logValues01,logValues02,logValues03)

submit <- combined_data[,c(1,2,3,4,17,18,19,20,33,34,35,36,5,6,7,8,21,22,23,24,37,38,39,40,9,10,11,12,25,26,27,28,41,42,43,44,13,14,15,16)]

colnames(submit) <- c("Sample_1","Sample_2","Sample_3","Sample_4","Sample_5","Sample_6","Sample_7","Sample_8","Sample_9","Sample_10","Sample_11","Sample_12","Sample_13","Sample_14","Sample_15","Sample_16","Sample_17","Sample_18","Sample_19","Sample_20","Sample_21","Sample_22","Sample_23","Sample_24","Sample_25","Sample_26","Sample_27","Sample_28","Sample_29","Sample_30","Sample_31","Sample_32","Sample_33","Sample_34","Sample_35","Sample_36","Sample_37","Sample_38","Sample_39","Sample_40")

write.table(submit,file="submit.txt")