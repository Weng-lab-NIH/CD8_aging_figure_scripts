
R_plotExploratory <- function(exploratory_plot, location, name) {
	print(paste('Plot was saved as an exploratory figure as ', name, sep = ''))
	print(paste('location: ', './figures/exploratory_figures/',location,'/',name,'.png',sep=''))
	png(filename=paste('./figures/exploratory_figures/',location,'/',name,'.png',sep=''), height = 1200, width = 1200, units='px')
	exploratory_plot
	dev.off()
}
