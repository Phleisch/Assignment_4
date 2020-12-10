import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("./result.csv")

# Plot
data = data.set_index("array_size")

data.plot(style=".")
plt.xlabel("Array size")
plt.ylabel("Execution time (s)")
plt.show()

# Format for latex table
pd.set_option('display.max_colwidth', 5000)
data["latex"] = data.index.to_series().apply(str) + " & " + data["cpu_time"].apply(str) + " & " + data["gpu_time"].apply(str) + "\\\\"
print(data["latex"].to_string(index=False))

