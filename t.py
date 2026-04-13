import polars as pl
import numpy as np
import matplotlib.pyplot as plt

df = pl.read_csv(r"C:\Users\Goob\Downloads\pastel_chart.csv")

df.columns = ['Participant', 'Age (years)', 'Sex', 'Age (months)', 
              'Time Lived With Both Eyes (months)', 'Mean Congruent Audiovisual Reaction Time', 
              'Mean Congruent Audiovisual Correct', 'Mean McGurk Reaction Time', 
              'McGurkBa', 'McGurkGa', 'McGurkDa']

dfBV = df.filter(pl.col("Participant").str.starts_with("BV"))
dfME = df.filter(pl.col("Participant").str.starts_with("ME"))

lineOfBestFitBV = np.polyfit(dfBV["Time Lived With Both Eyes (months)"], dfBV["McGurkDa"], 1)
lineOfBestFitME = np.polyfit(dfME["Time Lived With Both Eyes (months)"], dfME["McGurkDa"], 1)

plt.scatter(dfBV["Time Lived With Both Eyes (months)"], dfBV["McGurkDa"], label="BV", color="blue")
plt.plot(dfBV["Time Lived With Both Eyes (months)"], np.poly1d(lineOfBestFitBV)(dfBV["Time Lived With Both Eyes (months)"]), color="blue")
plt.scatter(dfME["Time Lived With Both Eyes (months)"], dfME["McGurkDa"], label="ME", color="red")
plt.plot(dfME["Time Lived With Both Eyes (months)"], np.poly1d(lineOfBestFitME)(dfME["Time Lived With Both Eyes (months)"]), color="red")
plt.xlabel("Time Lived With Both Eyes (months)")
plt.ylabel("McGurkDa")
plt.title("McGurkDa vs Time Lived With Both Eyes")
plt.legend()
plt.show()