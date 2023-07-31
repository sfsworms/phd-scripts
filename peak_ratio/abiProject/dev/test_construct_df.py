import pandas as pd


def my_function(dh, a=10):
    dh["A"] += a
    return dh


df = pd.DataFrame({"A":[7,8,9], "B":[9,10,11]})
print(df)

df = df.apply(my_function, a=20, axis=1)
print(df)
