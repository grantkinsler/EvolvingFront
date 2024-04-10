

# tools_path = f'{home_dir}/code/processing/tools/tools.py'
# sys.path.append(os.path.dirname(os.path.expanduser(tools_path)))
import tools as tools

from matplotlib import gridspec
from matplotlib.gridspec import GridSpec
import seaborn as sns
import matplotlib.pyplot as plt
# from matplotlib import colors
import matplotlib
import numpy as np
import pandas as p

import os

home_dir = '~/Documents/Stanford/Research/EvolvingFront/'
home_dir = os.path.expanduser(home_dir)


def plot_kdes(ax,this_anc,ancs,xdata,ydata,condense=False):

    if not condense:
        for anc in ancs:
            this_pure_diploid = this_anc[(this_anc['ancestor']==anc) & (this_anc['class_new']=='pure_diploids')]

            sns.kdeplot(x=this_pure_diploid[xdata].values,y=this_pure_diploid[ydata].values,
                        color=tools.anc_color_map[anc],alpha=0.4,thresh=0.2,levels=4,linestyles='--',ax=ax)

            this_neutral_haploid = this_anc[(this_anc['ancestor']==anc) & (this_anc['class_new']=='neutral_haploids')]

            sns.kdeplot(x=this_neutral_haploid[xdata].values,y=this_neutral_haploid[ydata].values,
                        color=tools.anc_color_map[anc],alpha=0.4,thresh=0.2,levels=4,ax=ax)
    else:
        this_pure_diploid = this_anc[(this_anc['ancestor'].isin(ancs)) & (this_anc['class_new']=='pure_diploids')]

        sns.kdeplot(x=this_pure_diploid[xdata].values,y=this_pure_diploid[ydata].values,
                    color='k',alpha=0.4,thresh=0.2,levels=4,linestyles='--',ax=ax)

        this_neutral_haploid = this_anc[(this_anc['ancestor'].isin(ancs)) & (this_anc['class_new']=='neutral_haploids')]

        sns.kdeplot(x=this_neutral_haploid[xdata].values,y=this_neutral_haploid[ydata].values,
                    color='k',alpha=0.4,thresh=0.2,levels=4,ax=ax)


    return ax



def tradeoff_figure(xdata,ydata,merged_fitness,
		ancestor_list=[['WT'],['CYR1','GPB2','TOR1','IRA1_MIS','IRA1_NON'],['WT'],['CYR1'],['GPB2'],['TOR1'],['IRA1_MIS'],['IRA1_NON']],
		evo_cond_list=['Evo2D','Evo3D'],
        pathway_list=['Ras/PKA','TOR/Sch9','HOG','RTG','TCA cycle','Deadenylation/Mitochondial Function'],
        # gene_exclusion_list=['MTH1'],
        centroid_cutoff = 2,
		centroids = False,
		pathways = True,
		annotate = False,
		nullcline_fitness = 'Fit2D_early_fitness',
		innovation_list = {},
        savefig=False,
        gray_alpha=0.2,
        legend=False,
        relative=False
		):


    fig = plt.figure(figsize=(12,14))
    outer_gs = gridspec.GridSpec(2, 1,height_ratios=[6,8])

    # larger gridspec for First-step and all second-step mutants
    few_inner_gs = gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec = outer_gs[0],wspace=0.25,hspace=0.25)
    # few_inner_gs = gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec = outer_gs[0],wspace=0,hspace=0)

    # smaller grid spec for each individual set of second-step mutants
    many_inner_gs = gridspec.GridSpecFromSubplotSpec(2,3,subplot_spec = outer_gs[1],wspace=0.25,hspace=0.25)

    for a,ancs in enumerate(ancestor_list):
        
        this_anc = merged_fitness[merged_fitness['ancestor'].isin(ancs)]
        this_anc = this_anc[this_anc['evolution_condition'].isin(evo_cond_list)]
        
        these_pure_diploids = this_anc[this_anc['class_new']=='pure_diploids']['barcode'].values
        these_neutral_haploids = this_anc[this_anc['class_new']=='neutral_haploids']['barcode'].values
        
        interesting_muts = this_anc[~this_anc['barcode'].isin(list(these_neutral_haploids)+list(these_pure_diploids))]
        
        if len(ancs) > 1:
            ax = fig.add_subplot(few_inner_gs[1]) 
            plot_kdes(ax,this_anc,ancs,xdata+'_relative',ydata+'_relative',condense=True)
            
        elif (ancs[0] == 'WT') & (a == 0):
            ax = fig.add_subplot(few_inner_gs[0]) 
            plot_kdes(ax,this_anc,ancs,xdata,ydata)
            
        else:
            ax = fig.add_subplot(many_inner_gs[a-2])
            plot_kdes(ax,this_anc,ancs,xdata,ydata)

#         for evo_cond in np.unique(this_anc['evolution_condition'].values):
        for ploidy,ploidy_list in {'haploid':['Haploid','haploid',np.nan,'?','NotSequenced','other'],'diploid':['diploid','Diploid']}.items():
#         for evo_cond in evo_cond_list:

            gene_list = {}
#             doubles_list = []

            this_data = interesting_muts[interesting_muts['ploidy_new'].isin(ploidy_list)]
            # gray_alpha = 0.1
            if centroids:
                bold_alpha = 0.3
            else:
                bold_alpha = 0.7
#             if len(ancs) < 2:
            if True:
                colors = []
                annotation_list = []
                for e,gene in enumerate(this_data['gene'].values):
                    color_assigned = matplotlib.colors.to_rgba('gray',gray_alpha)
                    already_assigned = False
                    if (not p.isnull(gene)) and (gene != 'NotSequenced'):
                        if gene in tools.mutation_color_map.keys():
                            if pathways:
                                color_assigned = matplotlib.colors.to_rgba(tools.find_pathway_color(gene),bold_alpha)
                                pathway = tools.find_pathway(gene)
                                
                                if pathway in gene_list.keys():
                                    gene_list[pathway].append(e)
                                else:
                                    gene_list[pathway] = [e]

                            else:
                                color_assigned = matplotlib.colors.to_rgba(tools.mutation_color_map[gene],bold_alpha)
                            
                                if gene in gene_list.keys():
                                    gene_list[gene].append(e)
                                else:
                                    gene_list[gene] = [e]
                        elif '+' in gene:
                            color_assigned = matplotlib.colors.to_rgba(tools.mutation_color_map['double_mutant'],bold_alpha)
                            if annotate:
	                            if len(ancs) < 2 and ancs[0] != 'WT':
	                                if this_data['barcode'].values[e] in innovation_list[f'{xdata}_{ydata}'][f'{anc}']:
	                                    annotation_list.append([e,gene])
#                                 else:
#                                     annotation_list.append([e,gene])
#                             if gene in gene_list.keys():
#                                 gene_list[gene].append(e)
#                             else:
#                                 gene_list[gene] = [e]
                        else:
                            color_assigned = matplotlib.colors.to_rgba('gray',gray_alpha)
                    else:
                            color_assigned = matplotlib.colors.to_rgba('gray',gray_alpha)
                    colors.append(color_assigned)
            else:
                colors = [tools.color_map[anc][evo_cond] for anc,evo_cond in zip(this_data['ancestor'],this_data['evolution_condition'])]

            if len(ancs) > 1:

                if centroids:
                    alpha = 0.7

                    plt.scatter(this_data[xdata+'_relative'].values,this_data[ydata+'_relative'].values,linewidths=0,alpha=gray_alpha,
                            color=colors,marker=tools.ploidy_marker_map[ploidy],s=15,label=f'{ancs[0]} {ploidy}')
                    for gene,e_list in gene_list.items():

                        if len(e_list) > centroid_cutoff:

                            gene_centroid = tools.centroid(this_data[[xdata+'_relative',ydata+'_relative']].values[e_list,:])

                            plt.scatter(gene_centroid[0],gene_centroid[1],
                                    color=colors[e_list[0]],edgecolors='k',linewidth=0.5,
                                        marker=tools.ploidy_marker_map[ploidy],s=40,alpha=0.9)

                else:
                    plt.scatter(this_data[xdata+'_relative'].values,this_data[ydata+'_relative'].values,linewidths=0,
                            color=colors,marker=tools.ploidy_marker_map[ploidy])

                plt.title('All Second Step')
               
                
            else:
                if centroids:
                    alpha = 0.7

                    plt.scatter(this_data[xdata].values,this_data[ydata].values,linewidths=0,alpha=gray_alpha,
                            color=colors,marker=tools.ploidy_marker_map[ploidy],s=15,label=f'{ancs[0]} {ploidy}')
                    for gene,e_list in gene_list.items():

                        # if len(e_list) > centroid_cutoff:

                        gene_centroid = tools.centroid(this_data[[xdata,ydata]].values[e_list,:])

                        plt.scatter(gene_centroid[0],gene_centroid[1],
                            color=colors[e_list[0]],edgecolors='k',linewidth=0.5,
                                marker=tools.ploidy_marker_map[ploidy],s=40,alpha=0.9)
                else:
                    alpha = 0.7

                    plt.scatter(this_data[xdata].values,this_data[ydata].values,linewidths=0,
                            color=colors,marker=tools.ploidy_marker_map[ploidy])

                for e,doubles in annotation_list:
                    if annotate:
                        plt.annotate(text=doubles.replace('+','+\n'),xy=(this_data[xdata].values[e],this_data[ydata].values[e]),
                                 xytext=(this_data[xdata].values[e]+0.001,this_data[ydata].values[e]+0.001),
                                 color='k',fontsize=8)
                    plt.scatter(this_data[xdata].values[e],this_data[ydata].values[e],
                                marker=tools.ploidy_marker_map[ploidy],
                                color='k',s=15,alpha=0.5,linewidth=0.2)

                # plt.legend(loc='lower left',fontsize='small')
                plt.title(f'{ancs[0]}')
                

        if len(ancs) < 2:
            for anc in ancs:
                if anc != 'WT':
                    background_mutant = merged_fitness[merged_fitness['barcode']==tools.rebarcoding_source_mutants[anc]]
                    plt.scatter(background_mutant[xdata].values,background_mutant[ydata].values,
                                    marker='+',color=tools.anc_color_map[anc],s=100)

                    plt.axvline(merged_fitness[merged_fitness['ancestor']==anc][xdata+'_ancestor'].values[0],color=tools.anc_color_map[anc],alpha=0.2,zorder=0)
                    plt.axhline(merged_fitness[merged_fitness['ancestor']==anc][ydata+'_ancestor'].values[0],color=tools.anc_color_map[anc],alpha=0.2,zorder=0)



        if relative:
            plt.xlim(tools.lims[xdata+'_relative'][0],tools.lims[xdata+'_relative'][1])
            plt.ylim(tools.lims[ydata+'_relative'][0],tools.lims[ydata+'_relative'][1])
        else:
            plt.xlim(tools.lims[xdata][0],tools.lims[xdata][1])
            plt.ylim(tools.lims[ydata][0],tools.lims[ydata][1])


        if len(ancs) >1:

        	plt.xlabel(f'{tools.labels[xdata]} relative to parental strain')
	        plt.ylabel(f'{tools.labels[ydata]} relative to parental strain')

	#         plt.axvline(0,color='k',linestyle=':')
	#         plt.axhline(0,color='k',linestyle=':')
	        
	        plt.axvline(0,color='gray',linestyle='-',zorder=0)
	        plt.axhline(0,color='gray',linestyle='-',zorder=0)

       	else:

            plt.xlabel(tools.labels[xdata])
            plt.ylabel(tools.labels[ydata])

            #         plt.axvline(0,color='k',linestyle=':')
            #         plt.axhline(0,color='k',linestyle=':')

            plt.axvline(0,color='gray',linestyle=':',zorder=0)
            plt.axhline(0,color='gray',linestyle=':',zorder=0)

            if (xdata == 'FerPerHour') & (ydata == 'ResPerHour'):
                if relative:
                    f_list = [1.0,1.5,2.0,2.5,3.0]
                    ferms = np.linspace(tools.lims[xdata+'_relative'][0],tools.lims[xdata+'_relative'][1],100)
                else:
                    f_list = [1.0,1.5,2.0,2.5,3.,3.50]
                    ferms = np.linspace(tools.lims[xdata][0],tools.lims[xdata][1],100)


                for fitness in f_list:

                # # for fitness in [1.0,1.5,2.0,2.5,3.0,3.5]:
                # #     # ferms = np.linspace(tools.lims[xdata][0],tools.lims[xdata][1],100)
                # for fitness in [1.0,1.5,2.0,2.5,3.0]:
                #     ferms = np.linspace(tools.lims[xdata+'_relative'][0],tools.lims[xdata+'_relative'][1],100)

                    resps = (fitness-16*ferms)/28 # 2day = 16*F
                    
                    norm = matplotlib.colors.Normalize(vmin=np.nanmin(merged_fitness[nullcline_fitness]),
                            vmax=np.nanmax(merged_fitness[nullcline_fitness]))

                    cm = matplotlib.cm.ScalarMappable(norm=norm, cmap='Reds') 
                    
                    plt.plot(ferms,resps,color=cm.to_rgba(fitness),alpha=0.05)
                    plt.text(x=ferms[-1],y=resps[-1],s=f'{fitness}',color=cm.to_rgba(fitness),alpha=0.1,ha='right',va='top')

        if legend: 
            if len(ancs) > 1:
                if pathways:
                    pathway_labels = []
                    for g,gene in enumerate(tools.mutation_color_map.keys()):
                        if gene in tools.pathway_gene_map.keys():
                            pathway_labels.append(tools.pathway_gene_map[gene])
                        else:
                            pathway_labels.append(gene)
                    
    #                 for p,pathway in enumerate(pathway_labels):
    #                     plt.text(y=0.4-0.02*(p%(len(pathway_labels))/2)),x=0.01+0.1*int(p/len(pathway_labels)*2),s=f'{gene}',color=tools.find_pathway_color[pathway],transform=plt.gca().transAxes)

                else:
                    for g,gene in enumerate(tools.mutation_color_map.keys()):
                        plt.text(y=0.4-0.02*(g%(len(tools.mutation_color_map.keys())/2)),x=0.01+0.1*int(g/len(tools.mutation_color_map.keys())*2),s=f'{gene}',color=tools.mutation_color_map[gene],transform=plt.gca().transAxes)

    plt.tight_layout()

    if savefig:
        if relative:
            plt.savefig(f'{home_dir}/figures/analysis/tradeoffs/tradeoffs_{xdata}_{ydata}_relative_{"pathway" if pathways else "mutation"}_colors{"_centroids" if centroids else ""}{"_annotatedInnovations" if annotate else ""}.pdf',bbox_inches='tight')
        else:
            plt.savefig(f'{home_dir}/figures/analysis/tradeoffs/tradeoffs_{xdata}_{ydata}_{"pathway" if pathways else "mutation"}_colors{"_centroids" if centroids else ""}{"_annotatedInnovations" if annotate else ""}.pdf',bbox_inches='tight')

    return fig

def tradeoff_figure_by_pathway(xdata,ydata,merged_fitness,
        ancestor_list=[['WT'],['CYR1','GPB2','TOR1','IRA1_MIS','IRA1_NON']],
        evo_cond_list=['Evo2D','Evo3D'],
        pathway_list=['Ras/PKA','TOR/Sch9','HOG','RTG','TCA cycle','Deadenylation/Mitochondial Function','Others'],
        centroids = False,
        pathways = True,
        annotate = False,
        nullcline_fitness = 'Fit2D_early_fitness',
        innovation_list = {},
        savefig=False,
        gray_alpha=0.1,
        bold_alpha=0.5,
        centroid_alpha=0.9,
        relative=False,
        ):

    # if relative:
    #     xdata = xdata + '_relative'
    #     ydata = ydata + '_relative'
    
    

    for a,ancs in enumerate(ancestor_list):
        for pathway in pathway_list:
            fig = plt.figure(figsize=(4,4))
            
            this_anc = merged_fitness[merged_fitness['ancestor'].isin(ancs)]
            this_anc = this_anc[this_anc['evolution_condition'].isin(evo_cond_list)]
            
            these_pure_diploids = this_anc[this_anc['class_new']=='pure_diploids']['barcode'].values
            these_neutral_haploids = this_anc[this_anc['class_new']=='neutral_haploids']['barcode'].values
            
            interesting_muts = this_anc[~this_anc['barcode'].isin(list(these_neutral_haploids)+list(these_pure_diploids))]
            
            # if len(ancs) > 1:
            if ancs[0] != 'WT':
                ax = plt.gca()
                # ax = fig.add_subplot(few_inner_gs[1]) 
                plot_kdes(ax,this_anc,ancs,xdata+'_relative',ydata+'_relative',condense=True)
                
            elif ancs[0] == 'WT':
                # ax = fig.add_subplot(few_inner_gs[0]) 
                ax = plt.gca()
                plot_kdes(ax,this_anc,ancs,xdata,ydata)
                
            # else:
            #     # ax = fig.add_subplot(many_inner_gs[a-1])
            #     ax = plt.gca()
            #     plot_kdes(ax,this_anc,ancs,xdata,ydata)

            for ploidy,ploidy_list in {'haploid':['Haploid','haploid',np.nan,'?','NotSequenced','other'],'diploid':['diploid','Diploid']}.items():

                gene_list = {}

                this_data = interesting_muts[interesting_muts['ploidy_new'].isin(ploidy_list)]
                gray_alpha = 0.1
                if centroids:
                    bold_alpha = bold_alpha
                else:
                    bold_alpha = bold_alpha
    #             if len(ancs) < 2:
                if True:
                    colors = []
                    annotation_list = []
                    for e,gene in enumerate(this_data['gene'].values):
                        color_assigned = matplotlib.colors.to_rgba('gray',0.1)
                        already_assigned = False
                        if (not p.isnull(gene)) and (gene != 'NotSequenced'):
                            # if gene in tools.mutation_color_map.keys():
                            if gene in tools.gene_pathway_map[pathway]:
                                if pathways:
                                    color_assigned = matplotlib.colors.to_rgba(tools.find_pathway_color(gene),bold_alpha)
                                    pathway = tools.find_pathway(gene)
                                    
                                    if pathway in gene_list.keys():
                                        gene_list[pathway].append(e)
                                    else:
                                        gene_list[pathway] = [e]

                                else:
                                    color_assigned = matplotlib.colors.to_rgba(tools.mutation_color_map[gene],bold_alpha)
                                
                                    if gene in gene_list.keys():
                                        gene_list[gene].append(e)
                                    else:
                                        gene_list[gene] = [e]
                            elif '+' in gene:
                                color_assigned = matplotlib.colors.to_rgba(tools.mutation_color_map['double_mutant'],gray_alpha)
                                if annotate:
                                    if len(ancs) < 2 and ancs[0] != 'WT':
                                        if this_data['barcode'].values[e] in innovation_list[f'{xdata}_{ydata}'][f'{anc}']:
                                            annotation_list.append([e,gene])
                            else:
                                color_assigned = matplotlib.colors.to_rgba('gray',gray_alpha)
                        else:
                                color_assigned = matplotlib.colors.to_rgba('gray',gray_alpha)
                        colors.append(color_assigned)
                else:
                    colors = [tools.color_map[anc][evo_cond] for anc,evo_cond in zip(this_data['ancestor'],this_data['evolution_condition'])]

                # if len(ancs) > 1:

                if centroids:
                    # alpha = bold_alpha

                    # if relative:

                    plt.scatter(this_data[xdata+'_relative'].values,this_data[ydata+'_relative'].values,linewidths=0,
                            color=colors,marker=tools.ploidy_marker_map[ploidy],s=15,label=f'{ancs[0]} {ploidy}')
                    
                    for gene,e_list in gene_list.items():

                        gene_centroid = tools.centroid(this_data[[xdata+'_relative',ydata+'_relative']].values[e_list,:])

                        plt.scatter(gene_centroid[0],gene_centroid[1],
                                color=colors[e_list[0]],edgecolors='k',linewidth=0.5,
                                    marker=tools.ploidy_marker_map[ploidy],s=40,alpha=centroid_alpha,zorder=1000)

                else:
                    plt.scatter(this_data[xdata+'_relative'].values,this_data[ydata+'_relative'].values,linewidths=0,
                            color=colors,marker=tools.ploidy_marker_map[ploidy])

                plt.title(f'{" ".join(ancs)}')
                   
                    
                # else:
                #     if centroids:
                #         alpha = 0.3

                #         plt.scatter(this_data[xdata].values,this_data[ydata].values,linewidths=0,
                #                 color=colors,marker=tools.ploidy_marker_map[ploidy],s=15,label=f'{ancs[0]} {ploidy}')
                #         for gene,e_list in gene_list.items():

                #             gene_centroid = tools.centroid(this_data[[xdata,ydata]].values[e_list,:])

                #             plt.scatter(gene_centroid[0],gene_centroid[1],
                #                     color=colors[e_list[0]],edgecolors='k',linewidth=0.5,
                #                         marker=tools.ploidy_marker_map[ploidy],s=40,alpha=0.9)
                #     else:
                #         alpha = 0.7

                #         plt.scatter(this_data[xdata].values,this_data[ydata].values,linewidths=0,
                #                 color=colors,marker=tools.ploidy_marker_map[ploidy])

        
                #     for e,doubles in annotation_list:
                #         if annotate:
                #             plt.annotate(text=doubles.replace('+','+\n'),xy=(this_data[xdata].values[e],this_data[ydata].values[e]),
                #                      xytext=(this_data[xdata].values[e]+0.001,this_data[ydata].values[e]+0.001),
                #                      color='k',fontsize=8)
                #         plt.scatter(this_data[xdata].values[e],this_data[ydata].values[e],
                #                     marker=tools.ploidy_marker_map[ploidy],
                #                     color='k',s=15,alpha=0.5,linewidth=0.2)

                #     # plt.legend(loc='lower left',fontsize='small')
                #     plt.title(f'{" ".join(ancs)}')
                    

            # if len(ancs) < 2:
            #     for anc in ancs:
            #         if anc != 'WT':
            #             background_mutant = merged_fitness[merged_fitness['barcode']==tools.rebarcoding_source_mutants[anc]]
            #             plt.scatter(background_mutant[xdata].values,background_mutant[ydata].values,
            #                             marker='+',color=tools.anc_color_map[anc],s=100)

            #             plt.axvline(merged_fitness[merged_fitness['ancestor']==anc][xdata+'_ancestor'].values[0],color=tools.anc_color_map[anc],alpha=0.2,zorder=0)
            #             plt.axhline(merged_fitness[merged_fitness['ancestor']==anc][ydata+'_ancestor'].values[0],color=tools.anc_color_map[anc],alpha=0.2,zorder=0)

            

            if relative:

                plt.xlim(tools.lims[xdata+'_relative'][0],tools.lims[xdata+'_relative'][1])
                plt.ylim(tools.lims[ydata+'_relative'][0],tools.lims[ydata+'_relative'][1])
            else:
                plt.xlim(tools.lims[xdata][0],tools.lims[xdata][1])
                plt.ylim(tools.lims[ydata][0],tools.lims[ydata][1])
            

            # if len(ancs) >1:
            if relative:

                plt.xlabel(f'{tools.labels[xdata]} relative to parental strain')
                plt.ylabel(f'{tools.labels[ydata]} relative to parental strain')

        #         plt.axvline(0,color='k',linestyle=':')
        #         plt.axhline(0,color='k',linestyle=':')
                
                plt.axvline(0,color='gray',linestyle='-',zorder=0)
                plt.axhline(0,color='gray',linestyle='-',zorder=0)

            else:

                plt.xlabel(tools.labels[xdata])
                plt.ylabel(tools.labels[ydata])

                #         plt.axvline(0,color='k',linestyle=':')
                #         plt.axhline(0,color='k',linestyle=':')

                plt.axvline(0,color='gray',linestyle=':',zorder=0)
                plt.axhline(0,color='gray',linestyle=':',zorder=0)

                if (xdata == 'FerPerHour') & (ydata == 'ResPerHour'):
                    if relative:
                        f_list = [1.0,1.5,2.0,2.5,3.0]
                        ferms = np.linspace(tools.lims[xdata+'_relative'][0],tools.lims[xdata+'_relative'][1],100)
                    else:
                        f_list = [1.0,1.5,2.0,2.5,3.,3.50]
                        ferms = np.linspace(tools.lims[xdata][0],tools.lims[xdata][1],100)


                    for fitness in f_list:
                    
                        resps = (fitness-16*ferms)/28 # 2day = 16*F
                        
                        norm = matplotlib.colors.Normalize(vmin=np.nanmin(merged_fitness[nullcline_fitness]),
                                vmax=np.nanmax(merged_fitness[nullcline_fitness]))

                        cm = matplotlib.cm.ScalarMappable(norm=norm, cmap='Reds') 
                        
                        plt.plot(ferms,resps,color=cm.to_rgba(fitness),alpha=0.05)
                        plt.text(x=ferms[-1],y=resps[-1],s=f'{fitness}',color=cm.to_rgba(fitness),alpha=0.1,ha='right',va='top')

            
            if len(ancs) > 1:
                if pathways:
                    pathway_labels = []
                    for g,gene in enumerate(tools.mutation_color_map.keys()):
                        if gene in tools.pathway_gene_map.keys():
                            pathway_labels.append(tools.pathway_gene_map[gene])
                        else:
                            pathway_labels.append(gene)
                    
    #                 for p,pathway in enumerate(pathway_labels):
    #                     plt.text(y=0.4-0.02*(p%(len(pathway_labels))/2)),x=0.01+0.1*int(p/len(pathway_labels)*2),s=f'{gene}',color=tools.find_pathway_color[pathway],transform=plt.gca().transAxes)

                # else:
                #     for g,gene in enumerate(tools.mutation_color_map.keys()):
                #         plt.text(y=0.4-0.02*(g%(len(tools.mutation_color_map.keys())/2)),x=0.01+0.1*int(g/len(tools.mutation_color_map.keys())*2),s=f'{gene}',color=tools.mutation_color_map[gene],transform=plt.gca().transAxes)

            plt.tight_layout()

            if savefig:
                if relative:
                    plt.savefig(f'{home_dir}/figures/analysis/tradeoffs/tradeoffs_{xdata}_{ydata}_relative_{"-".join(ancs)}_{pathway.replace("/","")}.pdf',bbox_inches='tight')
                else:
                    plt.savefig(f'{home_dir}/figures/analysis/tradeoffs/tradeoffs_{xdata}_{ydata}_{"-".join(ancs)}_{pathway.replace("/","")}.pdf',bbox_inches='tight')

    # return fig
