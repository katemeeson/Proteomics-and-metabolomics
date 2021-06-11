

def do_pca(X_df, n_components=20, centered=True):

    import numpy as np
    import pandas as pd
    from sklearn.decomposition import PCA
    
    '''
        Return:
            W: PCA loadings for each feature (how much each feature contributes to each PC)
            Xproj: projected data in PCA space (component scores, reduced dimension representation of data)
            fracs: Fraction of data explained by each principal component in a vector
    '''
    n_features = X_df.shape[1] # Original data dimension (no of features)
    if n_components > n_features:
        n_components = n_features
    if isinstance(X_df, pd.DataFrame):
        X = X_df.values # Get the data matrix from the X_df dataframe
    if centered:            
        X = X.astype('float')  # Since X is object
        X = X - X.mean(0)
        X = X/X.std(0)
    pca = PCA(n_components=n_components) # run the PCA algorithm from sklearn on data X
    pca.fit(X)
    fracs = pca.explained_variance_ratio_  # vector with explained variance for each principal component (PC)
    Xproj = pca.fit_transform(X)   # Low-dim projection (aka Scores) - n_sample, n_redim
    W_l = pca.components_            # PC Loadings - n_redim, n_feature
    # construct two DataFrames for later use
    # Loadings in a dataframe
    W_df = pd.DataFrame(W_l, index=['PC' + str(i) for i in np.arange(1, W_l.shape[0]+1)], columns=X_df.columns.values)
    # Low-dimensional projection of data in a dataframe
    Xproj_df = pd.DataFrame(Xproj, index=X_df.index, columns=['PC'+str(i) for i in np.arange(1, Xproj.shape[1]+1)])
    W = W_df.T
    return W, Xproj_df, fracs

def pca_biplot(W, scores, data, topN=None, XPC='PC1', YPC='PC2', feature=None):

	import plotly
	import plotly.graph_objs as go
	import plotly.express as px
	import pandas as pd
	import numpy as np
	
	#combined = pd.concat([scores, data], axis=1, sort=False)
    
	FeaturePCMagnitudes = W.loc[:,XPC]**2+W.loc[:,YPC]**2 # Squared length of vectors in the 2D PCA space
	FeaturePCMagnitudes.sort_values(inplace=True, ascending=False) # Sort from largest to smallest
	topN = 15 # Only show the 15 largest feature vectors

	if feature == None:
		fig = px.scatter(scores, x=XPC, y=YPC, color='label') # default is to plot the classes
		fig.update_traces(mode='markers', marker_line_width=1, marker_size=8)
	else:                               # otherwise plot the feature given in the function argument
		fig = go.Figure()
		tmp = np.array([feature + " "]*len(data))
		markertext = list()
		for i in range (1,len(data)):
			markertext.append(tmp[i] + str(data[feature].values[i]))
		fig.add_trace(go.Scatter(
			name="Data point",
    		x=scores.loc[:,XPC],
    		y=scores.loc[:,YPC],
    		text = markertext,
    		marker=dict(
    			line_width=1,
        		size=8,
         		color=data.loc[:,feature],
         		colorbar=dict(
             		title=feature,
             		x=-0.2,
             		y=0.5
         		),
         		colorscale="Viridis"
     		),
    		mode="markers"))
		fig.update_layout(xaxis_title=XPC, yaxis_title=YPC)

	for i, g in enumerate(FeaturePCMagnitudes.iteritems()): 
		if(topN is not None and i >= topN):
			continue
		fig.add_scatter(x=[0,W.loc[g[0],XPC]*2.5],y=[0,W.loc[g[0],YPC]*2.5],mode='lines',name=g[0])

	return fig

