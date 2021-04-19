 #!/bin/bash
 # 1. Go to https://disc.gsfc.nasa.gov/datasets?keywords=imerg&page=1
 # 2. Click Subset/Get data
 # 3. Click get data after applying filters and download the subset...txt file 
 wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --auth-no-challenge=on --keep-session-cookies --content-disposition -i subset_GPM_3IMERGHH_06_20210411_212529.txt
