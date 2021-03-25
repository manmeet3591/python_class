fig,ax = plt.subplots(ncols=3,nrows=3, figsize=(11.69,8.27), subplot_kw={'projection': ccrs.PlateCarree()})
cnt = 0
for i in range(3):
    for j in range(3):
        #print(pr_aer[0].sel(time=))
        tci_control = compute_tci(mrsos_control[cnt].sel(time=mrsos_control[cnt].time.dt.month.isin([5,6, 7, 8,9])), \
                                  hfls_control[cnt].sel(time=mrsos_control[cnt].time.dt.month.isin([5,6, 7, 8,9])))
        tci_aer = compute_tci(mrsos_aer[cnt].sel(time=mrsos_aer[cnt].time.dt.month.isin([5,6, 7, 8,9])), \
                                  hfls_aer[cnt].sel(time=mrsos_aer[cnt].time.dt.month.isin([5,6, 7, 8,9])))
#         tci_control = compute_tci(mrsos_control[cnt], hfls_control[cnt])
#         tci_aer = compute_tci(mrsos_aer[cnt], hfls_aer[cnt])
        
        diff_ = tci_aer - tci_control
        #diff_ = (pr_aer[cnt]-pr_control[cnt]).sel(time=pr_aer[cnt].time.dt.month.isin([5,6, 7, 8,9])).mean(dim='time')*86400
        diff_.sel(lon=slice(60,100),lat=slice(5,40)).plot(ax=ax[i,j], \
                                                                cmap='BrBG', extend='both', vmin=-10, vmax=10)
        ax[i,j].coastlines()
        ax[i,j].set_title(titles[cnt])
        cnt=cnt+1
        #print(i,j)
        if cnt>=6:
            break
        

fig.delaxes(ax[2,1])
fig.delaxes(ax[2,2])
fig.patch.set_linewidth(10)
fig.patch.set_edgecolor('cornflowerblue')
fig.suptitle('Terrestrial coupling index (Dirmeyer et al 2011) (AER - Control) MJJAS', size=16)
