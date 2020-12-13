array_sizes=( 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 120000, 140000, 160000, 180000, 200000, 250000, 300000, 350000, 400000, 450000, 500000, 550000, 600000)

outfile="result.csv"

echo "array_size,cpu_time,gpu_time" > ${outfile}

for size in "${array_sizes[@]}"
do
  ./exercise_2.out ${size} >> ${outfile}
done
