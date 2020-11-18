#!/usr/bin/env python3 


import json
from uk_covid19 import Cov19API
import pandas as pd 
import plotly.express as px

all_nations = [
    "areaType=nation"
]

cases_and_deaths = {
    "date": "date",
    "areaName": "areaName",
    "areaCode": "areaCode",
    "newCasesByPublishDate": "newCasesByPublishDate",
    "cumCasesByPublishDate": "cumCasesByPublishDate",
    "newDeathsByDeathDate": "newDeathsByDeathDate",
    "cumDeathsByDeathDate": "cumDeathsByDeathDate"
}


api_c_d = Cov19API(filters=all_nations, structure=cases_and_deaths)

data_c_d_json = api_c_d.get_json()

df = pd.read_json(json.dumps(data_c_d_json['data']) , orient='list')
df = df.loc[(df['date'] < '2020-01-03') | (df['date'] > '2020-03-01')]
data_wales = df[df['areaName'] == 'Wales']

fig_wales = px.bar(data_wales, x='date', y='newCasesByPublishDate', labels={'newCasesByPublishDate': 'Number of Cases', 'date': 'Date'})
fig_wales.update_traces(marker_color='#850404')
fig_wales.update_layout(title_text='Total number of SARS-CoV2 cases per day Wales')
fig_wales.show()


data_england = df[df['areaName'] == 'England']
fig_england = px.bar(data_england, x='date', y='newCasesByPublishDate', labels={'newCasesByPublishDate': 'Number of Cases', 'date': 'Date'})
fig_england.update_traces(marker_color='#017801')
fig_england.update_layout(title_text='Total number of SARS-CoV2 cases per day England')
fig_england.show()


data_scotland = df[df['areaName'] == 'Scotland']
fig_scotland = px.bar(data_scotland, x='date', y='newCasesByPublishDate', labels={'newCasesByPublishDate': 'Number of Cases', 'date': 'Date'})
fig_scotland.update_traces(marker_color='#013D78')
fig_scotland.update_layout(title_text='Total number of SARS-CoV2 cases per day Scotland')
fig_scotland.show()
