function(el, x) {
  // when hovering over an element, do something
  el.on('plotly_hover', function(d) {

    // extract tooltip text
    txt = d.points[0].data.text;
    // image is stored locally
    image_location = '/Users/pma/Dropbox/git_repos/mapp-metabolomics-unit/biostat_toolbox/data/' + txt + '.jpg';

    // define image to be shown
    var img = {
      // location of image
      source: image_location,
      // top-left corner
      x: 0,
      y: 1,
      sizex: 0.2,
      sizey: 0.2,
      xref: 'paper',
      yref: 'paper'
    };

    // show image and annotation 
    d3.select(el.id, {
        images: [img] 
    });
  })
}
