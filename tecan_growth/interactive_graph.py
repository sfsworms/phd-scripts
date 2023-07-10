# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 11:42:56 2023

@author: worms
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import matplotlib.pyplot as plt

def plot_growth_curve(file_path):
    # Load the spreadsheet
    xls = pd.ExcelFile(file_path)
    
    # Load the specified sheet into a dataframe
    df = xls.parse('Sheet0')

    # Identify the row containing the well names
    well_names_row_index = df[df.applymap(lambda x: x == 'A1' if pd.notnull(x) else False)].dropna(how='all').index[0]

    # Read the data again, but this time skip rows until the well names
    df_data = xls.parse('Sheet0', skiprows=well_names_row_index)

    # Drop unnecessary columns and set appropriate column names
    df_data.columns = df_data.iloc[0]
    df_data = df_data.drop(0)

    # Reset the index for the dataframe
    df_data = df_data.reset_index(drop=True)

    # Remove rows where "Cycle Nr." is not a number
    df_data = df_data[pd.to_numeric(df_data['Cycle Nr.'], errors='coerce').notna()]

    # Convert the "Cycle Nr.", "Time [s]" and "Temp. [°C]" columns to numeric
    df_data['Cycle Nr.'] = pd.to_numeric(df_data['Cycle Nr.'])
    df_data['Time [s]'] = pd.to_numeric(df_data['Time [s]'])
    df_data['Temp. [°C]'] = pd.to_numeric(df_data['Temp. [°C]'])

    # Convert all well data to numeric
    for col in df_data.columns[3:]:
        df_data[col] = pd.to_numeric(df_data[col])

    # Define the strains
    strains = ["A", "B", "D", "E", "F", "G", "H"]

    # Create a color palette from green to red
    colors = plt.cm.RdYlGn_r(np.linspace(0, 1, len(strains)))

    # Create a new plot
    fig = go.Figure()

    # Add a trace for each strain
    for strain, color in zip(strains, colors):
        # Select columns for the current strain
        strain_cols = [col for col in df_data.columns if col.startswith(strain)]
        # Calculate the average OD and standard deviation
        avg_OD = df_data[strain_cols].mean(axis=1)
        std_OD = df_data[strain_cols].std(axis=1)
        # Add the average OD trace
        fig.add_trace(go.Scatter(
            x=df_data['Time [s]'], y=avg_OD,
            name=strain,
            line=dict(color='rgba'+str(tuple(int(c*255) for c in color[:3]))),
            mode='lines'
        ))
        # Add the standard deviation trace
        fig.add_trace(go.Scatter(
            x=df_data['Time [s]'], y=avg_OD + std_OD,
            name=strain+' upper',
            mode='lines',
            marker=dict(color="#444"),
            line=dict(width=0),
            showlegend=False
        ))
        fig.add_trace(go.Scatter(
            x=df_data['Time [s]'], y=avg_OD - std_OD,
            name=strain+' lower',
            marker=dict(color="#444"),
            line=dict(width=0),
            mode='lines',
            fillcolor='rgba'+str(tuple(int(c*255) for c in color[:3])+(0.1,)),
            fill='tonexty',
            showlegend=False
        ))

    # Set plot title and labels
    fig.update_layout(
        title='Average OD for each strain over time',
        xaxis_title='Time [s]',
        yaxis_title='Average OD',
    )

    # Show the plot
    fig.show()


# Use the function
plot_growth_curve(r"C:\Users\worms\Dropbox\data\tecan\2023.06.22_test_select_no_ara\2023.06.22 test growth tecan.xlsx")
