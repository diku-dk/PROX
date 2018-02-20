from ProxXML import RawData, load_from_xml


if __name__ == '__main__':

    data = RawData()

    load_from_xml('/Users/kenny/Documents/PROX_PAPERS/2016_CHUNKED_BVH/2017_RESULTS/raw_data/funnel_w8100/scene.xml', data, verbose=True)
    load_from_xml('/Users/kenny/Documents/PROX_PAPERS/2016_CHUNKED_BVH/2017_RESULTS/raw_data/glasses_w8100/scene.xml', data, verbose=True)
    #load_from_xml('/Users/kenny/Documents/PROX_PAPERS/2016_CHUNKED_BVH/2017_RESULTS/raw_data/temple_w8100/scene.xml', data, verbose=True)
    #load_from_xml('/Users/kenny/Documents/PROX_PAPERS/2016_CHUNKED_BVH/2017_RESULTS/raw_data/colosseum_w8100/scene.xml', data, verbose=True)
