import matplotlib.pyplot as plt
import pandas as pd


class Sample:
    def __init__(self, name, x, y):
        self.name = name
        self.x = x
        self.y = y


def plot(df):
    plt.plot(df["READ"].y, df["READ"].x)
    plt.title(df["READ"].name)
    plt.savefig(df["READ"].name)
    plt.close()


datas = []
datas.append(Sample("read1", [0,1,2,3,4,5,6], [4,5,6,9,8,7,4]))
datas.append(Sample("read2", [2,1,2,3,4,5,6], [4,5,6,9,8,7,4]))
datas.append(Sample("read3", [4,6,2,3,4,5,6], [4,5,6,9,8,7,4]))

df = pd.DataFrame()
for sample in datas:
    df = pd.concat([df, pd.DataFrame({"READ": [sample]})], ignore_index=True)

df.apply(plot, 1)
