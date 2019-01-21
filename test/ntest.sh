export Z3_FMM=0
export Z4_GEO=0
export Z5_SMT=0
export Z6_EQU=0
rm stop *.o *.out
make clean
touch output
while [ $Z3_FMM -le 3 ]
do
echo "||-----------------------------------------------------------------||"
echo "||   fmm :  $Z3_FMM     geo :  $Z4_GEO     smt :  $Z5_SMT     equ :  $Z6_EQU      compile  ||"
echo "||-----------------------------------------------------------------||"
make ex_n
make par
export Z6_EQU=`expr $Z6_EQU + 1`
if [ $Z6_EQU -gt 0 ]; then
export Z6_EQU=0
export Z5_SMT=`expr $Z5_SMT + 1`
fi
if [ $Z5_SMT -gt 1 ]; then
export Z5_SMT=0
export Z4_GEO=`expr $Z4_GEO + 1`
fi
if [ $Z4_GEO -gt 3 ]; then
export Z4_GEO=0
export Z3_FMM=`expr $Z3_FMM + 1`
fi
if [ -e "stop" ]; then
echo "||-----------------------------------------------------------------||"
echo "||                        canceled by user                         ||"
echo "||-----------------------------------------------------------------||"
export Z3_FMM=4
fi
done
rm *.o *.out
make clean
