from sklearn.cluster import KMeans, DBSCAN, AffinityPropagation
from sklearn import metrics


def PCA_process(X, nps):
    from sklearn.decomposition import PCA
    print('Shape of data to PCA:', X.shape)
    pca = PCA(n_components=nps)
    X_PC = pca.fit_transform(X)     
    print('Shape of data output by PCA:', X_PC.shape)
    print('PCA recover:', pca.explained_variance_ratio_.sum())
    return X_PC



def Kmeans_cluster(X_embedding, n_clusters, merge=False):
    cluster_model = KMeans(n_clusters=n_clusters, init='k-means++', n_init=100, max_iter=1000, tol=1e-6)
    cluster_labels = cluster_model.fit_predict(X_embedding)
    if merge:
        cluster_labels = merge_cluser(X_embedding, cluster_labels)
    score = metrics.silhouette_score(X_embedding, cluster_labels, metric='euclidean')
    
    return cluster_labels, score

def merge_cluser(X_embedding, cluster_labels):
    count_dict, out_count_dict = {}, {}
    for cluster in cluster_labels:
        count_dict[cluster] = count_dict.get(cluster, 0) + 1
    clusters = count_dict.keys()
    n_clusters = len(clusters)
    for cluster in clusters: 
        out_count_dict[cluster] = count_dict[cluster] 
    for cluster in clusters: 
        cur_n = count_dict[cluster]
        if cur_n <=3:
            min_dis = 1000
            merge_to = cluster
            center_cluster = X_embedding[cluster_labels==cluster].mean(0)
            for cluster_2 in clusters:
                if cluster_2 == cluster:
                    continue
                center_cluster_2 = X_embedding[cluster_labels==cluster_2].mean(0)
                dist = np.linalg.norm(center_cluster - center_cluster_2)
                if dist < min_dis:
                    min_dis = dist
                    merge_to = cluster_2

            cluster_labels[cluster_labels==cluster] = merge_to
            print('Merge group', cluster, 'to group', merge_to, 'with', cur_n, 'samples')
            out_count_dict[cluster] = 0
            out_count_dict[merge_to] += cur_n
            if cluster < n_clusters-1:
                cluster_labels[cluster_labels==n_clusters-1] = cluster
                print('Group', n_clusters-1, 'is renamed to group', cluster)
                out_count_dict[cluster] = out_count_dict[n_clusters-1]
                del out_count_dict[n_clusters-1]
            print(out_count_dict)

    return cluster_labels