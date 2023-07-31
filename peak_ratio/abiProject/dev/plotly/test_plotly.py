from Bio import SeqIO
import matplotlib.pyplot as plt
import plotly.express as px


def plot(data_list, filename):
    fig = px.line(y=data_list, title=filename)
    fig.write_html('first_figure.html', auto_open=True)


def extract_raw_data(filename):
    rec = SeqIO.read(filename, "abi")
    return [rec.annotations["abif_raw"]["DATA1"], rec.annotations["abif_raw"]["DATA2"], rec.annotations["abif_raw"]["DATA3"], rec.annotations["abif_raw"]["DATA4"]]


def run():
    while True:
        filename = input("Type filename here: ")
        data = extract_raw_data(filename)
        plot(data, filename)


if __name__ == "__main__":
    run()
