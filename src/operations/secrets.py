"""
This file was created by the CMU DELPHI Epi-forecasting group.
https://github.com/cmu-delphi/

Contains various non-public strings.
"""


class db:
    auto = ('{SECRET_DB_USERNAME_AUTO}', '{SECRET_DB_PASSWORD_AUTO}')
    backup = ('{SECRET_DB_USERNAME_BACKUP}', '{SECRET_DB_PASSWORD_BACKUP}')
    epi = ('{SECRET_DB_USERNAME_EPI}', '{SECRET_DB_PASSWORD_EPI}')


class api:
    twitter = '{SECRET_API_AUTH_TWITTER}'
    ght = '{SECRET_API_AUTH_GHT}'
    fluview = '{SECRET_API_AUTH_FLUVIEW}'
    cdc = '{SECRET_API_AUTH_CDC}'
    quidel = '{SECRET_API_AUTH_QUIDEL}'
    norostat = '{SECRET_API_AUTH_NOROSTAT}'
    afhsb = '{SECRET_API_AUTH_AFHSB}'
    sensors = '{SECRET_API_AUTH_SENSORS}'

    class sensor_subsets:
        twtr_sensor = '{SECRET_API_AUTH_TWTR_SENSOR}'
        gft_sensor = '{SECRET_API_AUTH_GFT_SENSOR}'
        ght_sensors = '{SECRET_API_AUTH_GHT_SENSORS}'
        cdc_sensor = '{SECRET_API_AUTH_CDC_SENSOR}'
        quid_sensor = '{SECRET_API_AUTH_QUID_SENSOR}'
        wiki_sensor = '{SECRET_API_AUTH_WIKI_SENSOR}'
