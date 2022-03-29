import cdsapi

c = cdsapi.Client()
for year in range(1979,2020):
  for i_mon in range(1,13,3):
    c.retrieve(
        'sis-agrometeorological-indicators',
        {
            'format': 'tgz',
            'variable': 'solar_radiation_flux',
            'year': str(year),
            'month': [
                '0'+str(i_mon), '0'+str(i_mon+1), '0'+str(i_mon+2),
            ],
            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31',
            ],
        },
        'agera5_'+str(year)+'_'+str(i_mon)+'_'+str(i_mon+1)+'.tar.gz')
