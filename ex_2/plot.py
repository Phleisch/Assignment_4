import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("./result.csv")

# Plot
data = data.set_index("array_size")
data["cpu_time"] = data["cpu_time"] * 1000
data["gpu_time"] = data["gpu_time"] * 1000

data.plot(style=".")
plt.xlabel("Array size")
plt.ylabel("Execution time (ms)")
plt.show()

data["gpu_time"].plot(style=".")
plt.xlabel("Array size")
plt.ylabel("GPU Execution time (ms)")
plt.show()

# Format for latex table
data["cpu_time"] = data["cpu_time"].round(4)
data["gpu_time"] = data["gpu_time"].round(4)
pd.set_option('display.max_colwidth', 5000)
data["latex"] = data.index.to_series().apply(str) + " & " + data["cpu_time"].apply(str) + " & " + data["gpu_time"].apply(str) + "\\\\"
print(data["latex"].to_string(index=False))

