/* Browser-only explorer: saved views read published outputs; Advanced downloads plans without modifying batch data. */
(() => {
  'use strict';
  const $ = (id) => document.getElementById(id);
  const tabs = ['overview', 'pca', 'pcoa', 'volcano', 'plsda', 'features', 'classes', 'downloads', 'advanced'];
  const hex = /^[0-9a-f]{64}$/;
  const outputPath = /^(tables\/[A-Za-z0-9_.-]+\.tsv|plots\/[A-Za-z0-9_.-]+\.(?:png|pdf))$/;
  const palette = ['#197b8a', '#ca7544', '#664ba0', '#42805f', '#a94564', '#976522', '#3459a2', '#77853b'];
  let catalogue = [];
  let builderSource = null;
  let current = null;
  let generation = 0;
  const detachedMode = new URL(location.href).searchParams.get('detached') === '1';
  let activeTab = detachedMode ? 'features' : 'overview';
  let metadata = [];
  let sampleGroupColumns = [];
  let defaultGroup = '';
  let selectedFeature = '';
  let featureCandidate = 0;

  let groupColors = new Map();
  let drawerVersion = 0;
  let detachedWindow = null;

  function el(tag, attrs = {}, ...children) {
    const node = document.createElement(tag);
    for (const [key, value] of Object.entries(attrs)) {
      if (key === 'class') node.className = value;
      else if (key === 'text') node.textContent = value;
      else if (key === 'onClick') node.addEventListener('click', value);
      else node.setAttribute(key, value);
    }
    node.append(...children.filter((child) => child !== null && child !== undefined));
    return node;
  }
  const text = (value) => String(value ?? '');
  const label = (value) => value === '' || value === null || value === undefined || value === 'NA' ? 'Not recorded' : text(value);
  const number = (value) => value === '' || value === 'NA' ? NaN : Number(value);
  const escapePlot = (value) => text(value).replace(/[&<>"']/g, (char) => ({'&':'&amp;', '<':'&lt;', '>':'&gt;', '"':'&quot;', "'":'&#39;'}[char]));
  function message(container, value, error = false) {
    container.replaceChildren(el('p', {class: error ? 'notice error' : 'notice', text: value}));
  }
  function status(value, error = false) {
    $('status').textContent = value;
    $('status').classList.toggle('error', error);
  }
  function valid(ctx) { return ctx === current && ctx.version === generation; }
  function panel(name) { return $(name + '-content'); }
  function card(title, ...children) { return el('section', {class:'card'}, el('h3', {text:title}), ...children); }
  function note(value, error = false) { return el('p', {class: error ? 'notice error' : 'note', text:value}); }
  function control(title, id, values, initial) {
    const select = el('select', {id});
    for (const [value, name] of values) select.append(el('option', {value, text:name}));
    if (initial !== undefined) select.value = initial;
    return el('div', {}, el('label', {for:id, text:title}), select);
  }
  function appendLinkedText(cell, value) {
    const content = label(value);
    const urls = /https?:\/\/[^\s<>"';]+/gi;
    let start = 0, match, linked = false;
    while ((match = urls.exec(content))) {
      const candidate = match[0].replace(/[.,)]+$/, '');
      let url;
      try { url = new URL(candidate); } catch { continue; }
      if (!['http:', 'https:'].includes(url.protocol)) continue;
      cell.append(document.createTextNode(content.slice(start,match.index)),
        el('a',{href:url.href,target:'_blank',rel:'noopener noreferrer',referrerpolicy:'no-referrer',text:candidate}));
      start = match.index + candidate.length;
      linked = true;
    }
    if (linked) {
      cell.append(document.createTextNode(content.slice(start)));
      cell.classList.add('linked-value');
    } else cell.textContent = content;
  }
  function table(headers, rows, onFeature, filename = 'dashboard-table.tsv') {
    const size = 100;
    const visible = new Set(headers.map((_,i) => i));
    const filters = headers.map(() => '');
    const numericColumns = headers.map((_,i) => {
      let found = false;
      for (const row of rows) {
        const value = row[i];
        if (value === null || value === undefined || value === '' || value === 'NA') continue;
        if (!Number.isFinite(number(value))) return false;
        found = true;
      }
      return found;
    });
    const numericLiteral = '(-?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:e[+-]?\\d+)?)';
    const threshold = new RegExp('^([<>]=?|=)\\s*'+numericLiteral+'$','i');
    const interval = new RegExp('^'+numericLiteral+'\\s*\\.\\.\\s*'+numericLiteral+'$','i');
    function matcher(query, numeric) {
      if (!query) return null;
      if (numeric) {
        const comparison = threshold.exec(query.trim());
        if (comparison) {
          const limit = Number(comparison[2]);
          const operator = comparison[1];
          return (value) => {
            if (value === null || value === undefined) return false;
            const n = number(value);
            return Number.isFinite(n) && (operator === '>' ? n > limit :
              operator === '>=' ? n >= limit : operator === '<' ? n < limit :
                operator === '<=' ? n <= limit : n === limit);
          };
        }
        const bounds = interval.exec(query.trim());
        if (bounds) {
          const lower = Number(bounds[1]), upper = Number(bounds[2]);
          return (value) => {
            if (value === null || value === undefined) return false;
            const n = number(value);
            return Number.isFinite(n) && n >= lower && n <= upper;
          };
        }
      }
      const needle = query.toLocaleLowerCase();
      return (value) => label(value).toLocaleLowerCase().includes(needle);
    }
    const collator = new Intl.Collator(undefined,{numeric:true,sensitivity:'base'});
    let sortIndex = -1, ascending = true, page = 0, displayed = [];
    const body = el('tbody');
    const columnOptions = el('div',{class:'table-column-options'});
    const columns = el('details',{class:'table-columns'},el('summary',{text:'Columns'}),columnOptions);
    const checkboxes = headers.map((heading,i) => {
      const checkbox = el('input',{type:'checkbox'});
      checkbox.checked = true;
      checkbox.addEventListener('change',() => {
        if (checkbox.checked) visible.add(i);
        else visible.delete(i);
        page = 0;
        update();
      });
      columnOptions.append(el('label',{},checkbox,document.createTextNode(heading)));
      return checkbox;
    });
    const downloadButton = el('button',{type:'button',class:'button secondary small',text:'Download displayed TSV'});
    downloadButton.addEventListener('click',() => {
      const indexes = headers.map((_,i) => i).filter((i) => visible.has(i));
      const quote = (value) => {
        const field = text(value);
        return /[\t\r\n"]/.test(field) ? '"' + field.replaceAll('"','""') + '"' : field;
      };
      const lines = [indexes.map((i) => quote(headers[i])).join('\t'),
        ...displayed.map((row) => indexes.map((i) => quote(label(row[i]))).join('\t'))];
      const url = URL.createObjectURL(new Blob([lines.join('\n')+'\n'],{type:'text/tab-separated-values;charset=utf-8'}));
      const link = el('a',{href:url,download:filename});
      document.body.append(link);
      link.click();
      link.remove();
      setTimeout(() => URL.revokeObjectURL(url),60000);
    });
    const headings = headers.map((heading,i) => {
      const sort = el('button',{type:'button',class:'table-sort',text:heading+' ↕',
        'aria-label':'Sort by '+heading});
      sort.addEventListener('click',() => {
        ascending = sortIndex === i ? !ascending : true;
        sortIndex = i;
        page = 0;
        update();
      });
      const filter = el('input',{type:'search',
        placeholder:numericColumns[i] ? 'Filter… >2 or 1..3' : 'Filter…',
        title:numericColumns[i] ? 'Numeric: >2, >=2, <2, <=2, =2, or inclusive 1..3; otherwise match text' : 'Match text',
        'aria-label':'Filter '+heading});
      filter.addEventListener('input',() => { filters[i] = filter.value; page = 0; update(); });
      return el('th',{scope:'col'},sort,filter);
    });
    const wrap = el('div',{class:'table-wrap'},el('table',{},el('thead',{},el('tr',{},...headings)),body));
    const previous = el('button',{type:'button',class:'button secondary small',text:'Previous',onClick:() => { page--; render(); }});
    const next = el('button',{type:'button',class:'button secondary small',text:'Next',onClick:() => { page++; render(); }});
    const count = el('span');
    const pager = el('div',{class:'pager'},previous,count,next);
    function render() {
      const end = Math.min((page+1)*size,displayed.length);
      const fragment = document.createDocumentFragment();
      for (const row of displayed.slice(page*size,end)) {
        const tr = el('tr');
        headers.forEach((heading,i) => {
          if (!visible.has(i)) return;
          const cell = el('td',{title:label(row[i])});
          if (onFeature && i === 0 && row[i])
            cell.append(el('button',{type:'button',class:'text-button',text:row[i],onClick:() => onFeature(row[i])}));
          else appendLinkedText(cell,row[i]);
          tr.append(cell);
        });
        fragment.append(tr);
      }
      body.replaceChildren(fragment);
      previous.disabled = page === 0;
      next.disabled = end === displayed.length;
      count.textContent = displayed.length ? (page*size+1)+'–'+end+' / '+displayed.length+' rows' : '0 / 0 rows';
    }
    function update() {
      for (let i = 0; i < headers.length; i++) {
        headings[i].hidden = !visible.has(i);
        headings[i].setAttribute('aria-sort',sortIndex === i ? (ascending ? 'ascending' : 'descending') : 'none');
        headings[i].firstChild.textContent = headers[i] + (sortIndex === i ? (ascending ? ' ↑' : ' ↓') : ' ↕');
        checkboxes[i].disabled = visible.size === 1 && visible.has(i) || Boolean(onFeature && i === 0);
      }
      const matchers = filters.map((query,i) => matcher(query,numericColumns[i]));
      displayed = rows.filter((row) => matchers.every((match,i) => !match || match(row[i])));
      if (sortIndex !== -1) displayed.sort((a,b) => {
        const left = label(a[sortIndex]), right = label(b[sortIndex]);
        const x = number(a[sortIndex]), y = number(b[sortIndex]);
        const result = Number.isFinite(x) && Number.isFinite(y) ? x-y : collator.compare(left,right);
        return ascending ? result : -result;
      });
      render();
    }
    update();
    return el('div',{class:'data-table'},
      el('div',{class:'table-tools'},columns,downloadButton),wrap,pager);
  }
  function metrics(entries) {
    return el('div', {class:'metrics'}, ...entries.map(([name, value]) => el('div', {class:'metric'}, el('strong', {text:label(value)}), el('span', {text:name}))));
  }

  // A TSV field may contain tabs, escaped quotes, or quoted CR/LF. Empty fields remain empty strings.
  function parseTSV(source) {
    const records = [];
    let record = [], field = '', quoted = false, closed = false;
    for (let i = 0; i < source.length; i++) {
      const c = source[i];
      if (quoted) {
        if (c === '"' && source[i + 1] === '"') { field += '"'; i++; }
        else if (c === '"') { quoted = false; closed = true; }
        else field += c;
      } else if (c === '"' && field === '' && !closed) quoted = true;
      else if (c === '\t') { record.push(field); field = ''; closed = false; }
      else if (c === '\n' || c === '\r') {
        record.push(field); records.push(record);
        record = []; field = ''; closed = false;
        if (c === '\r' && source[i + 1] === '\n') i++;
      } else if (closed || c === '"') throw new Error('Malformed quoted TSV field at character ' + i);
      else field += c;
    }
    if (quoted) throw new Error('Unterminated quoted TSV field');
    if (field !== '' || record.length || closed) records.push([...record, field]);
    if (!records.length) throw new Error('Empty TSV file');
    const headers = records.shift();
    if (headers.length) headers[0] = headers[0].replace(/^\uFEFF/, '');
    return records.map((values, i) => {
      if (values.length !== headers.length) throw new Error('TSV column count differs at row ' + (i + 2));
      return Object.fromEntries(headers.map((header, j) => [header, values[j]]));
    });
  }
  function listed(ctx, path) {
    if (!outputPath.test(path) || !Object.hasOwn(ctx.run.outputs, path) || !hex.test(ctx.run.outputs[path]))
      throw new Error('File is not a valid published output: ' + path);
  }
  async function verifiedBytes(ctx, path) {
    listed(ctx, path);
    if (!window.crypto?.subtle) throw new Error('SHA-256 verification requires HTTPS or loopback HTTP (localhost). This page cannot use unverified files.');
    const response = await fetch('../' + ctx.run.hash + '/' + path, {signal:ctx.controller.signal});
    if (!response.ok) throw new Error('Could not load ' + path + ' (HTTP ' + response.status + ')');
    const bytes = await response.arrayBuffer();
    const digest = Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256', bytes)), (x) => x.toString(16).padStart(2, '0')).join('');
    if (digest !== ctx.run.outputs[path]) throw new Error('SHA-256 mismatch: ' + path + '. The published file may be corrupt.');
    return bytes;
  }
  function saved(ctx, path, optional = false) {
    if (optional && !Object.hasOwn(ctx.run.outputs, path)) return Promise.resolve(null);
    if (!ctx.cache.has(path)) {
      const request = verifiedBytes(ctx, path).then((bytes) => parseTSV(new TextDecoder('utf-8', {fatal:true}).decode(bytes)));
      ctx.cache.set(path, request);
      request.catch(() => { if (ctx.cache.get(path) === request) ctx.cache.delete(path); });
    }
    return ctx.cache.get(path);
  }
  const rows = (ctx, name, optional = false) => saved(ctx, 'tables/' + name + '.tsv', optional);
  const preferredStructureSources = ['gnps','isdb','sirius'];
  function recordedText(value) {
    const string = text(value).trim();
    return ['','na','null','none'].includes(string.toLowerCase()) ? '' : string;
  }
  function structureSource(source, library = '') {
    const name = text(source).trim().toLowerCase();
    const savedLibrary = recordedText(library).toLowerCase();
    if (name === 'met_annot_enhancer') return savedLibrary === 'isdb' ? 'isdb' : '';
    if (name === 'isdb' && savedLibrary && savedLibrary !== 'isdb') return '';
    return preferredStructureSources.includes(name) ? name : '';
  }
  function structureKey(value) {
    const match = /^([a-z]{14})(?:-|$)/i.exec(text(value).trim());
    return match ? match[1].toUpperCase() : '';
  }
  function galleryStructures(candidates, horizontal) {
    const structures = new Map();
    function add(smiles, source, key, record) {
      const value = text(smiles).trim();
      if (!value || ['na','null','none'].includes(value.toLowerCase())) return;
      const ik2d = structureKey(key);
      // A shared SMILES is not evidence of a shared identity when the saved keys are missing.
      const group = source && ik2d ? 'key:'+ik2d :
        'unkeyed:'+record.provenance+'\u0000'+source+'\u0000'+value;
      if (!structures.has(group)) structures.set(group,{ik2d:source ? ik2d : '',records:[]});
      structures.get(group).records.push({...record,smiles:value,source,ik2d});
    }
    for (const row of candidates) {
      const source = structureSource(row.source,row.library);
      const provenance = source ? source.toUpperCase()+
        (source === 'isdb' && !recordedText(row.library) ? ' (library not recorded)' : '') :
        (recordedText(row.library) || recordedText(row.source) || 'Source not recorded');
      add(row.smiles,source,row.inchikey,{label:row.label,provenance,
        libraryKnown:!!recordedText(row.library),score:row.taxonomic_score,taxon:row.organism});
    }
    for (const row of horizontal) for (const sourceName of preferredStructureSources) {
      const library = sourceName === 'isdb' ? row.isdb_libname : '';
      const source = sourceName === 'isdb' && Object.hasOwn(row,'isdb_libname') && !recordedText(library) ?
        '' : structureSource(sourceName,library);
      const provenance = source ? source.toUpperCase()+
        (source === 'isdb' && !recordedText(library) ? ' (library not recorded)' : '') :
        (recordedText(library) || 'Source not recorded')+' (horizontal '+sourceName+'_SMILES)';
      add(row[sourceName+'_SMILES'],source,row[sourceName+'_IK2D'],
        {label:'Saved horizontal structure',provenance,
          libraryKnown:sourceName !== 'isdb' || !!recordedText(library),
          score:sourceName === 'isdb' ? row.isdb_score_taxo : '',
          taxon:sourceName === 'isdb' ? row.isdb_organism_name : ''});
    }
    return structures;
  }
  function sourceStar(source, records) {
    const related = records.filter((record) => record.source === source);
    const scored = source === 'isdb' && related.find((record) =>
      recordedText(record.score) && Number.isFinite(number(record.score)));
    const taxon = recordedText((scored && recordedText(scored.taxon) ? scored :
      related.find((record) => recordedText(record.taxon)))?.taxon);
    const score = scored ? Number(scored.score) : NaN;
    const description = source.toUpperCase()+' saved 2D-key source'+
      (source === 'isdb' ? '; taxonomic score '+
        (Number.isFinite(score) ? scored.score+' (visual color scale 0–8)' : 'not recorded')+
        '; reported taxon: '+(taxon || 'not recorded')+
        (related.some((record) => !record.libraryKnown) ? '; library provenance not recorded' : '')+
        '; score is not biological confirmation' : '')+'.';
    const star = el('span',{class:'source-star'+(source === 'isdb' ? ' source-star-isdb' : ' source-star-'+source),
      role:'img',title:description,'aria-label':description},el('span',{'aria-hidden':'true',text:'★'}),
      document.createTextNode(' '+source.toUpperCase()+
        (Number.isFinite(score) ? ' · '+scored.score : '')));
    if (Number.isFinite(score)) {
      const fraction = Math.min(8,Math.max(0,score))/8;
      star.style.backgroundColor = 'hsl('+(49+79*fraction)+' 81% '+(82-32*fraction)+'%)';
    }
    return star;
  }
  const smilesEncoder = new TextEncoder();
  const smilesBuffer = new Uint8Array(5001);
  let hoverFields = new Set(['row m/z','row retention time']);
  let hoverSource = '';
  let hoverCard = null;
  function depictionUrl(smiles, width, height) {
    return 'https://api.naturalproducts.net/latest/depict/2D?smiles='+
      encodeURIComponent(smiles)+'&toolkit=rdkit&width='+width+'&height='+height;
  }
  function usableSmiles(value) {
    const smiles = text(value).trim();
    if (!smiles || smiles.length > 5000 || ['na','null','none'].includes(smiles.toLowerCase())) return '';
    const result = smilesEncoder.encodeInto(smiles,smilesBuffer);
    return result.read === smiles.length && result.written <= 5000 ? smiles : '';
  }
  function hideFeatureHover() { hoverCard?.remove(); hoverCard = null; }
  function fillHoverSettings(index) {
    const fields = $('hover-field-choices');
    fields.replaceChildren(el('legend',{text:'Hover fields'}),
      el('label',{},el('input',{type:'checkbox',checked:'',disabled:''}),document.createTextNode('Feature ID (always shown)')));
    for (const name of index.fields) {
      const checkbox = el('input',{type:'checkbox',value:name});
      checkbox.checked = hoverFields.has(name);
      checkbox.addEventListener('change',() => {
        if (checkbox.checked) hoverFields.add(name);
        else hoverFields.delete(name);
        hideFeatureHover();
      });
      fields.append(el('label',{},checkbox,document.createTextNode(name)));
    }
    const source = $('hover-structure-source');
    source.replaceChildren(el('option',{value:'',text:'Any saved source'}),
      ...index.sources.map((name) => el('option',{value:name,text:name.toUpperCase()})));
    if (!index.sources.includes(hoverSource)) hoverSource = '';
    source.value = hoverSource;
    source.disabled = !index.sources.length;
  }
  function featureHoverIndex(ctx) {
    if (!ctx.hoverIndex) ctx.hoverIndex = Promise.all([
      rows(ctx,'variable_metadata',true),rows(ctx,'annotation_candidates',true),rows(ctx,'annotation_horizontal',true)
    ]).then(([variables,candidates,horizontal]) => {
      const features = new Map(), sources = new Set();
      const entry = (id) => {
        if (!features.has(id)) features.set(id,{fields:null,structures:new Map()});
        return features.get(id);
      };
      const addStructure = (id,source,value,library) => {
        const canonical = structureSource(source,library);
        if (!id || !canonical) return;
        const structures = entry(id).structures;
        if (structures.has(canonical)) return;
        const smiles = text(value).trim();
        if (smiles && !['na','null','none'].includes(smiles.toLowerCase())) {
          structures.set(canonical,smiles);
          sources.add(canonical);
        }
      };
      for (const row of variables ?? []) if (row.feature_id) entry(row.feature_id).fields = row;
      for (const row of horizontal ?? []) for (const source of preferredStructureSources)
        addStructure(row.feature_id,source,row[source+'_SMILES'],source === 'isdb' ? row.isdb_libname : '');
      for (const row of candidates ?? []) addStructure(row.feature_id,row.source,row.smiles,row.library);
      const fields = variables?.length ? Object.keys(variables[0]).filter((name) =>
        name && name !== 'feature_id' && variables.some((row) => row[name] && row[name] !== 'NA')) : [];
      const orderedSources = preferredStructureSources.filter((source) => sources.has(source));
      const index = {features,fields,sources:orderedSources};
      if (valid(ctx)) fillHoverSettings(index);
      return index;
    });
    return ctx.hoverIndex;
  }
  function bindFeatureHover(ctx,plot,index,featureId) {
    if (plot.featureHoverHandler && plot.removeListener) {
      plot.removeListener('plotly_hover',plot.featureHoverHandler);
      plot.removeListener('plotly_unhover',hideFeatureHover);
    }
    const onHover = (event) => {
      if (!valid(ctx) || !event.points?.length) return;
      const id = featureId(event.points[0]);
      if (!id) return;
      hideFeatureHover();
      const feature = index.features.get(id);
      const details = el('dl');
      for (const name of index.fields) if (hoverFields.has(name))
        details.append(el('dt',{text:name}),el('dd',{text:label(feature?.fields?.[name])}));
      const card = el('div',{class:'feature-hover',role:'tooltip'},el('strong',{text:'Feature '+id}),details);
      const source = hoverSource || index.sources.find((name) => feature?.structures.has(name));
      const smiles = feature?.structures.get(source);
      if (smiles) {
        card.append(el('p',{class:'note',text:'Saved '+source.toUpperCase()+' SMILES · candidate structure, not confirmed identity'}));
        if (usableSmiles(smiles)) {
          const image = el('img',{src:depictionUrl(smiles,300,200),referrerpolicy:'no-referrer',
            alt:'2D depiction of saved '+source.toUpperCase()+' structure for feature '+id});
          const fallback = note('Depiction unavailable from the external service.');
          fallback.hidden = true;
          image.addEventListener('error',() => { image.remove(); fallback.hidden = false; });
          card.append(image,fallback);
        } else card.append(note('Saved SMILES exceeds the depiction service’s 5000-byte limit.'));
      } else card.append(note('No saved '+(hoverSource ? hoverSource.toUpperCase() : 'structural')+' SMILES for this feature.'));
      document.body.append(card);
      hoverCard = card;
      const bounds = card.getBoundingClientRect();
      const x = event.event?.clientX ?? plot.getBoundingClientRect().left;
      const y = event.event?.clientY ?? plot.getBoundingClientRect().top;
      card.style.left = Math.max(8,Math.min(innerWidth-bounds.width-8,x+16))+'px';
      card.style.top = Math.max(8,Math.min(innerHeight-bounds.height-8,y+16))+'px';
    };
    plot.on('plotly_hover',onHover);
    plot.on('plotly_unhover',hideFeatureHover);
    plot.addEventListener('mouseleave',hideFeatureHover);
    plot.featureHoverHandler = onHover;
  }
  function featureColorData(ctx) {
    if (!ctx.featureColorData) ctx.featureColorData = Promise.all([
      rows(ctx,'variable_metadata',true),rows(ctx,'annotation_candidates',true),
      rows(ctx,'annotation_horizontal',true)
    ]).then(([variables,candidates,horizontal]) => {
      const candidateFields = [
        ['npc_pathway','CANOPUS NPC pathway'],['npc_superclass','CANOPUS NPC superclass'],
        ['npc_class','CANOPUS NPC class'],['source','Candidate record source (not 2D-key agreement)']
      ].filter(([key]) => candidates?.some((row) => row[key] && row[key] !== 'NA'))
        .map(([key,name]) => ({id:'candidate:'+key,key,name,rows:candidates}));
      const consensusFields = horizontal?.length ? [{
        id:'horizontal:sources_number_IK2D',key:'sources_number_IK2D',
        name:'Saved 2D-key agreement (0–3 sources)',rows:horizontal,
        numeric:true,range:[0,3],scaleLabel:'Agreeing sources',tickvals:[0,1,2,3]
      }] : [];
      const scoreFields = horizontal?.some((row) => Number.isFinite(number(row.isdb_score_taxo))) ? [{
        id:'horizontal:isdb_score_taxo',key:'isdb_score_taxo',
        name:'Saved ISDB taxonomic score (0–8)',rows:horizontal,
        numeric:true,range:[0,8],scaleLabel:'ISDB taxonomic score'
      }] : [];
      const measuredFields = [];
      for (const key of Object.keys(variables?.[0] ?? {})) {
        if (!key || key === 'feature_id' || key === 'row ID') continue;
        const values = [...new Set(variables.map((row) => row[key]).filter((value) => value && value !== 'NA'))];
        if (values.length < 2) continue;
        let numeric = true, min = Infinity, max = -Infinity;
        for (const value of values) {
          const parsed = number(value);
          if (!Number.isFinite(parsed)) { numeric = false; break; }
          min = Math.min(min,parsed);
          max = Math.max(max,parsed);
        }
        if (numeric && /\b(?:id|number)\b/i.test(key)) continue;
        if (numeric || values.length <= 60) measuredFields.push({
          id:'measured:'+key,key,name:'Measured · '+key,rows:variables,
          ...(numeric ? {numeric:true,range:[min,max]} : {})
        });
      }
      return {fields:[...consensusFields,...scoreFields,...candidateFields,...measuredFields],maps:new Map(),
        featureIds:(variables ?? []).map((row) => row.feature_id).filter(Boolean)};
    });
    return ctx.featureColorData;
  }
  function featureCategories(data,field) {
    if (!data.maps.has(field.id)) {
      const categories = new Map();
      for (const row of field.rows) {
        const value = text(row[field.key]).trim();
        if (!row.feature_id || !value || value === 'NA') continue;
        if (!categories.has(row.feature_id)) categories.set(row.feature_id,new Set());
        categories.get(row.feature_id).add(value);
      }
      data.maps.set(field.id,categories);
    }
    return data.maps.get(field.id);
  }
  function featureDisplayOptions(ctx,data,scope,ids,plot,onChange) {
    const choices = [['','Saved default'],...data.fields.map((field) => [field.id,field.name])];
    const fieldControl = control('Colour features by',scope+'-feature-colour',choices);
    const select = fieldControl.querySelector('select');
    const pickerList = el('div',{class:'group-color-choices'});
    const details = el('details',{class:'intensity-options'},el('summary',{text:'Display options · feature colours'}),
      el('div',{class:'intensity-options-body'},fieldControl,
        note('Numeric fields use a continuous Viridis scale; missing values stay “Not recorded” and conflicting saved values use “Multiple values”. Saved 2D-key counts come from horizontal evidence, not candidate records. Candidate record source lists SIRIUS, CANOPUS classes and met_annot_enhancer (ISDB/MS1 evidence), not GNPS agreement.'),pickerList));
    ctx.featureColorFields ??= new Map();
    ctx.featureCategoryColors ??= new Map();
    ctx.featureColorOpen ??= new Map();
    select.value = ctx.featureColorFields.get(scope) || '';
    details.open = ctx.featureColorOpen.get(scope) || false;
    details.addEventListener('toggle',() => ctx.featureColorOpen.set(scope,details.open));
    const field = () => data.fields.find((item) => item.id === select.value);
    const category = (id) => {
      const selected = field();
      if (!selected) return '';
      const values = featureCategories(data,selected).get(id);
      return !values?.size ? 'Not recorded' : values.size > 1 ? 'Multiple values' : values.values().next().value;
    };
    const color = (name) => ctx.featureCategoryColors.get(select.value+'\u0000'+name);
    function groups(rows, defaults) {
      const selected = field();
      if (!selected) return defaults();
      const categories = new Map();
      for (const row of rows) {
        const name = category(row.feature_id);
        if (!categories.has(name)) categories.set(name,[]);
        categories.get(name).push(row);
      }
      if (!selected.numeric) return [...categories].sort(([a],[b]) => a.localeCompare(b))
        .map(([name,subset]) => [name,subset,color(name),selected.id+'\u0000'+name]);
      const values = [];
      const other = [];
      for (const [name,subset] of categories) {
        if (Number.isFinite(number(name))) for (const row of subset) values.push(row);
        else other.push([name,subset,color(name),selected.id+'\u0000'+name]);
      }
      const scale = {colorscale:'Viridis',cmin:selected.range[0],cmax:selected.range[1],
        showscale:true,colorbar:{title:{text:escapePlot(selected.scaleLabel ?? selected.name)},
          ...(selected.tickvals ? {tickmode:'array',tickvals:selected.tickvals} : {})}};
      return [...(values.length ? [[selected.name,values,
        values.map((row) => number(category(row.feature_id))),'',scale,false]] : []),...other];
    }
    function updatePickers() {
      pickerList.replaceChildren();
      const selected = field();
      if (!selected) return;
      const categories = [...new Set(ids.map(category).filter((name) =>
        !selected.numeric || !Number.isFinite(number(name))))].sort((a,b) =>
        a === 'Not recorded' ? 1 : b === 'Not recorded' ? -1 : a.localeCompare(b));
      categories.forEach((name,i) => {
        const key = select.value+'\u0000'+name;
        if (!ctx.featureCategoryColors.has(key)) ctx.featureCategoryColors.set(key,
          name === 'Not recorded' ? '#a4a9af' : name === 'Multiple values' ? '#967d9d' : palette[i % palette.length]);
        const picker = el('input',{type:'color','aria-label':'Colour for '+name});
        picker.value = ctx.featureCategoryColors.get(key);
        picker.addEventListener('input',() => {
          ctx.featureCategoryColors.set(key,picker.value);
          plot.data?.forEach((trace,index) => {
            if (trace.legendgroup === key) Plotly.restyle(plot,{'marker.color':picker.value},[index]);
          });
        });
        pickerList.append(el('label',{},document.createTextNode(name),picker));
      });
    }
    select.addEventListener('change',() => {
      ctx.featureColorFields.set(scope,select.value);
      updatePickers();
      onChange();
    });
    updatePickers();
    return {node:details,field,category,color,groups};
  }

  function safeRender(name, renderer) {
    const ctx = current, target = panel(name);
    if (!ctx) return;
    const token = Symbol(name);
    ctx.viewToken = token;
    message(target, 'Loading verified saved results…');
    Promise.resolve().then(() => renderer(ctx, target, token)).catch((error) => {
      if (viewValid(ctx,name,token)) message(target, error.message || text(error), true);
    });
  }
  function viewValid(ctx,name,token) { return valid(ctx) && activeTab === name && ctx.viewToken === token; }
  function available(ctx, name) { return Object.hasOwn(ctx.run.outputs, 'tables/' + name + '.tsv'); }
  function savedPlots(ctx, stem) {
    return ['png','pdf'].map((extension) => 'plots/' + stem + '.' + extension)
      .filter((path) => Object.hasOwn(ctx.run.outputs,path));
  }
  function resizeActivePlots() {
    requestAnimationFrame(() => {
      if (window.Plotly && current) panel(activeTab).querySelectorAll('.plot.js-plotly-plot')
        .forEach((plot) => Plotly.Plots.resize(plot));
    });
  }
  function closeFeatureDrawer() {
    const drawer = $('feature-drawer');
    if (drawer.hidden) return;
    drawerVersion++;
    drawer.open = false;
    drawer.hidden = true;
    $('dashboard-main').classList.remove('drawer-open');
    if (window.Plotly) $('feature-drawer-content').querySelectorAll('.plot.js-plotly-plot')
      .forEach((plot) => Plotly.purge(plot));
    $('feature-drawer-content').replaceChildren();
    resizeActivePlots();
  }
  function syncDetached() {
    if (!detachedWindow) return false;
    if (detachedWindow.closed) {
      detachedWindow = null;
      $('detached-status').hidden = true;
      return false;
    }
    if (current) detachedWindow.postMessage({kind:'mapp-feature-sync',hash:current.run.hash,
      feature:selectedFeature,colors:[...groupColors]},location.origin);
    return true;
  }
  function publishFeatureSelection() {
    if (detachedMode) {
      $('detached-feature').value = selectedFeature;
      const url = new URL(location.href);
      if (selectedFeature) url.searchParams.set('feature',selectedFeature);
      else url.searchParams.delete('feature');
      history.replaceState(null,'',url);
      window.opener?.postMessage({kind:'mapp-feature-sync',hash:current.run.hash,feature:selectedFeature},location.origin);
    } else syncDetached();
  }
  function reattachFeature() {
    if (detachedWindow && !detachedWindow.closed) detachedWindow.close();
    detachedWindow = null;
    $('detached-status').hidden = true;
    if (selectedFeature && current && activeTab !== 'features') openFeature(selectedFeature);
    else if (activeTab === 'features') render();
  }
  function detachFeature() {
    if (!current || !selectedFeature || detachedMode) return;
    const url = new URL(location.href);
    url.searchParams.set('detached','1');
    url.searchParams.set('run',current.run.hash);
    url.searchParams.set('feature',selectedFeature);
    const popup = window.open(url.href,'_blank','popup,width=1120,height=900');
    if (!popup) {
      status('Popup blocked. Allow a new window for this dashboard, then choose Detach again.',true);
      return;
    }
    detachedWindow = popup;
    $('detached-status').hidden = false;
    closeFeatureDrawer();
    popup.focus();
  }

  function selectTab(name) {
    if (!tabs.includes(name)) return;
    if (name !== activeTab) closeFeatureDrawer();
    hideFeatureHover();
    activeTab = name;
    for (const tab of tabs) {
      $(tab).hidden = tab !== name;
      const button = document.querySelector('[data-tab="' + tab + '"]');
      button.classList.toggle('active', tab === name);
      if (tab === name) button.setAttribute('aria-current', 'page'); else button.removeAttribute('aria-current');
    }
    render();
  }
  function render() {
    if (!current) return;
    const dispatch = {overview:renderOverview, pca:(ctx,target,token) => renderOrdination(ctx,target,token,'pca','PC'),
      pcoa:(ctx,target,token) => renderOrdination(ctx,target,token,'pcoa','PCoA'), volcano:renderVolcano,
      plsda:renderPlsda, features:(ctx,target,token) => renderFeatures(ctx,target,token,detachedMode),
      classes:renderClasses, downloads:renderDownloads, advanced:renderAdvanced};
    safeRender(activeTab, dispatch[activeTab]);
  }
  function featurePlotKey(hash) { return 'mapp:feature-plot:' + hash; }
  function savedFeaturePlotSettings(ctx, stages, facets) {
    let saved;
    try { saved = JSON.parse(localStorage.getItem(featurePlotKey(ctx.run.hash))); } catch {}
    if (!saved || typeof saved !== 'object' || Array.isArray(saved)) saved = {};
    return {
      stage:stages.includes(saved.stage) ? saved.stage : stages.includes('transformed') ? 'transformed' : stages.includes('raw') ? 'raw' : stages[0],
      group:sampleGroupColumns.includes(saved.group) ? saved.group : defaultGroup,
      role:saved.role === '' || metadata.some((row) => row.sample_type === saved.role) ? saved.role : '',
      facet:facets.includes(saved.facet) ? saved.facet : '',
      groups:Array.isArray(saved.groups) && saved.groups.every((name) => typeof name === 'string') ? new Set(saved.groups) : null,
      open:saved.open === true
    };
  }
  function saveFeaturePlotSettings(ctx) {
    const settings = ctx.featurePlotSettings;
    try {
      const key = featurePlotKey(ctx.run.hash);
      const value = JSON.stringify({...settings,groups:settings.groups === null ? null : [...settings.groups]});
      if (localStorage.getItem(key) !== value) localStorage.setItem(key,value);
    } catch {}
  }
  function selector(select, values, placeholder) {
    select.replaceChildren();
    if (placeholder) select.append(el('option', {value:'', text:placeholder}));
    values.forEach((value) => select.append(el('option', {value, text:label(value)})));
  }
  function sampleKey(row) { return row.filename || row.sample || row.sample_id; }
  function groupName(row,field = defaultGroup) { return label(row[field]); }
  function groupColor(name) {
    if (!groupColors.has(name)) groupColors.set(name,palette[groupColors.size % palette.length]);
    return groupColors.get(name);
  }
  function updateGroupColor(name, color) {
    groupColors.set(name,color);
    for (const container of [panel(activeTab),$('feature-drawer-content')]) container.querySelectorAll('.plot[data-sample-groups]').forEach((plot) => {
      plot.data?.forEach((trace,index) => {
        if (trace.legendgroup !== name) return;
        Plotly.restyle(plot,trace.type === 'box' ?
          {'marker.color':color,'line.color':color} : {'marker.color':color},[index]);
      });
    });
    if (!detachedMode) syncDetached();
  }
  function setupMetadata(data) {
    metadata = data;
    if (!data.length) throw new Error('sample_metadata.tsv contains no samples');
    const forbidden = /^(sample|sample_id|filename)$|path|file|dir|uri|url/i;
    sampleGroupColumns = Object.keys(data[0]).filter((key) =>
      !forbidden.test(key) && new Set(data.map((row) => row[key]).filter(Boolean)).size <= 30);
    if (!sampleGroupColumns.length) throw new Error('No descriptive grouping columns in sample_metadata.tsv');
    defaultGroup = sampleGroupColumns.includes(current.run.design_group) ? current.run.design_group :
      sampleGroupColumns.includes('sample_type') ? 'sample_type' : sampleGroupColumns[0];
    groupColors = new Map();
    const primary = [...new Set(metadata.map((row) => groupName(row,defaultGroup)))].sort();
    const allNames = [...new Set(metadata.flatMap((row) => sampleGroupColumns.map((field) => groupName(row,field))))].sort();
    const choices = [...primary,...allNames.filter((name) => !primary.includes(name))];
    $('group-color-choices').replaceChildren(...choices.map((name) => {
      const picker = el('input',{type:'color','aria-label':'Colour for '+name});
      picker.value = groupColor(name);
      picker.addEventListener('input',() => updateGroupColor(name,picker.value));
      return el('label',{},document.createTextNode(name),picker);
    }));
  }
  function infoFacts(run) {
    const facts = [['Dataset',run.dataset_id], ['Run hash',run.hash], ['Description',run.description ?? 'Not in this catalogue'],
      ['Recipe SHA-256',run.recipe_sha256 ?? 'Not in this catalogue'],
      ['Changes vs previous dataset run',run.changes_vs_previous_dataset_run ?? 'Not in this catalogue'],
      ['Created',run.created_at], ['Saved schema',run.schema_version],
      ['Method','V2 saved analyses'], ['Design group',run.design_group], ['Intensity',run.intensity_measure], ['Transformation',run.transformation]];
    $('run-facts').replaceChildren(...facts.map(([key, value]) => el('div', {}, el('dt', {text:key}), el('dd', {text:label(value), class:key === 'Run hash' ? 'mono' : ''}))));
  }
  async function chooseRun(hash, updateHistory = false, feature = '') {
    if (updateHistory) {
      const url = new URL(location.href);
      url.searchParams.set('run', hash);
      history.pushState(null, '', url);
    }
    if (current) current.controller.abort();
    closeFeatureDrawer();
    hideFeatureHover();
    hoverFields = new Set(['row m/z','row retention time']);
    hoverSource = '';
    $('hover-field-choices').replaceChildren(el('legend',{text:'Hover fields'}),note('Open Settings or a feature plot to load saved fields.'));
    selector($('hover-structure-source'),[],'Any saved source');
    $('hover-structure-source').disabled = true;
    generation++;
    const run = catalogue.find((item) => item.hash === hash);
    current = null;
    metadata = [];
    selectedFeature = detachedMode ? feature : '';
    if (detachedMode) $('detached-feature').value = selectedFeature;
    featureCandidate = 0;
    $('group-color-choices').replaceChildren();
    for (const tab of tabs) {
      if (window.Plotly) panel(tab).querySelectorAll('.plot').forEach((node) => Plotly.purge(node));
      panel(tab).replaceChildren();
    }
    $('run-facts').replaceChildren();
    if (!run) {
      $('run-select').value = '';
      status('Unknown run hash in URL: ' + hash + '. Select an available run.', true);
      message(panel(activeTab), 'This run is not in the published catalogue. Choose an available run above.', true);
      return;
    }
    const ctx = {run, cache:new Map(), controller:new AbortController(), version:generation};
    current = ctx;
    if (!detachedMode) syncDetached();
    $('run-select').value = hash;
    infoFacts(run);
    status('Verifying sample metadata for ' + run.dataset_id + ' · ' + hash.slice(0,12) + '…');
    try {
      const samples = await rows(ctx, 'sample_metadata');
      if (!valid(ctx)) return;
      setupMetadata(samples);
      if (!detachedMode) syncDetached();
      status('Viewing ' + run.dataset_id + ' / ' + hash.slice(0,12) + ' · saved results only');
      render();
      if ($('hover-settings').open) featureHoverIndex(ctx).catch((error) => {
        if (valid(ctx)) message($('hover-field-choices'),error.message,true);
      });
    } catch (error) {
      if (valid(ctx)) { status(error.message, true); message(panel(activeTab), error.message, true); }
    }
  }
  async function boot() {
    try {
      const response = await fetch('runs.json', {cache:'no-store'});
      if (!response.ok) throw new Error('Could not load runs.json (HTTP ' + response.status + ')');
      const data = await response.json();
      if (data.schema_version !== 1 || !Array.isArray(data.runs)) throw new Error('Unsupported runs.json catalogue schema');
      builderSource = data.builder_source ?? null;
      catalogue = data.runs;
      if (!catalogue.length) throw new Error('No completed runs are published in runs.json');
      for (const run of catalogue) {
        if (!hex.test(run.hash) || !run.outputs || typeof run.outputs !== 'object' || Array.isArray(run.outputs) ||
            (run.description !== undefined && (typeof run.description !== 'string' || !run.description ||
              run.description.trim() !== run.description || [...run.description].length > 160 ||
              /[\x00-\x1f\x7f-\x9f]/u.test(run.description))) ||
            (run.recipe_sha256 !== undefined && (typeof run.recipe_sha256 !== 'string' || !hex.test(run.recipe_sha256))) ||
            (run.changes_vs_previous_dataset_run !== undefined &&
              (typeof run.changes_vs_previous_dataset_run !== 'string' ||
                !run.changes_vs_previous_dataset_run ||
                run.changes_vs_previous_dataset_run.trim() !== run.changes_vs_previous_dataset_run ||
                run.changes_vs_previous_dataset_run.length > 8192 ||
                /[\x00-\x1f\x7f-\x9f]/u.test(run.changes_vs_previous_dataset_run))) ||
            Object.entries(run.outputs).some(([path, sha]) => !outputPath.test(path) && path.startsWith('tables/') || !hex.test(sha))) {
          throw new Error('Invalid run entry or output checksum in runs.json');
        }
      }
      const select = $('run-select');
      select.replaceChildren(...catalogue.map((run, i) => el('option', {value:run.hash,
        text:run.hash.slice(0,12) + ' · ' + (run.description ?? run.dataset_id) +
          ' · ' + label(run.created_at) + (i === 0 ? ' · latest' : '')})));
      select.disabled = false;
      const requested = new URL(location.href).searchParams.get('run');
      await chooseRun(requested === null ? catalogue[0].hash : requested,false,
        detachedMode ? new URL(location.href).searchParams.get('feature') || '' : '');
      if (detachedMode) window.opener?.postMessage({kind:'mapp-feature-ready'},location.origin);
    } catch (error) {
      status(error.message, true);
      message(panel(activeTab), error.message, true);
    }
  }
  $('hover-settings').addEventListener('toggle',() => {
    const ctx = current;
    if ($('hover-settings').open && ctx) featureHoverIndex(ctx).catch((error) => {
      if (valid(ctx)) message($('hover-field-choices'),error.message,true);
    });
  });
  $('hover-structure-source').addEventListener('change',(event) => {
    hoverSource = event.target.value;
    hideFeatureHover();
  });
  $('run-select').addEventListener('change', (event) => chooseRun(event.target.value, true));
  window.addEventListener('popstate', () => {
    const url = new URL(location.href);
    chooseRun(url.searchParams.get('run') ?? catalogue[0]?.hash,false,
      detachedMode ? url.searchParams.get('feature') || '' : '');
  });
  window.addEventListener('storage',(event) => {
    if (!current || event.key !== featurePlotKey(current.run.hash)) return;
    current.featurePlotSettings = undefined;
    if (!selectedFeature) return;
    if (activeTab === 'features') render();
    else if ($('feature-drawer').open && !$('feature-drawer').hidden) {
      const ctx = current, version = ++drawerVersion;
      renderFeatures(ctx,$('feature-drawer-content'),ctx.viewToken,true,version).catch((error) => {
        if (valid(ctx) && version === drawerVersion) message($('feature-drawer-content'),error.message || text(error),true);
      });
    }
  });
  document.querySelectorAll('[data-tab]').forEach((button) => button.addEventListener('click', () => selectTab(button.dataset.tab)));
  $('feature-drawer').addEventListener('toggle',() => {
    const open = $('feature-drawer').open && !$('feature-drawer').hidden;
    $('dashboard-main').classList.toggle('drawer-open',open);
    resizeActivePlots();
  });
  $('feature-full').addEventListener('click',() => selectTab('features'));
  $('feature-detach').addEventListener('click',detachFeature);
  $('focus-detached').addEventListener('click',() => {
    if (!syncDetached()) reattachFeature();
    else detachedWindow.focus();
  });
  $('reattach-feature').addEventListener('click',reattachFeature);
  if (!detachedMode) window.addEventListener('focus',() => {
    if (detachedWindow?.closed) syncDetached();
  });
  $('detached-attach').addEventListener('click',() => {
    window.opener?.postMessage({kind:'mapp-feature-attach'},location.origin);
    window.close();
  });
  $('detached-feature').addEventListener('change',() => {
    const id = $('detached-feature').value.trim();
    if (id) openFeature(id);
    else if (detachedMode && current) {
      selectedFeature = '';
      publishFeatureSelection();
      render();
    }
  });
  window.addEventListener('message',(event) => {
    if (event.origin !== location.origin || !event.data || typeof event.data !== 'object') return;
    const data = event.data;
    if (detachedMode) {
      if (!window.opener || event.source !== window.opener || data.kind !== 'mapp-feature-sync' ||
          !hex.test(data.hash) || typeof data.feature !== 'string' || data.feature.length > 256) return;
      const url = new URL(location.href);
      url.searchParams.set('run',data.hash);
      if (data.feature) url.searchParams.set('feature',data.feature);
      else url.searchParams.delete('feature');
      history.replaceState(null,'',url);
      const applyColors = () => {
        for (const entry of Array.isArray(data.colors) ? data.colors : []) {
          if (Array.isArray(entry) && typeof entry[0] === 'string' && entry[0].length <= 128 &&
              typeof entry[1] === 'string' && /^#[0-9a-fA-F]{6}$/.test(entry[1]) &&
              groupColors.get(entry[0]) !== entry[1]) updateGroupColor(entry[0],entry[1]);
        }
      };
      if (current?.run.hash !== data.hash) {
        if (!catalogue.some((run) => run.hash === data.hash)) {
          status('Selected run is absent from this window’s catalogue. Republish and reload the detached window.',true);
          return;
        }
        chooseRun(data.hash,false,data.feature).then(() => {
          if (current?.run.hash === data.hash) applyColors();
        });
      } else {
        if (selectedFeature !== data.feature) {
          selectedFeature = data.feature;
          featureCandidate = 0;
          $('detached-feature').value = selectedFeature;
          if (metadata.length) render();
        }
        applyColors();
      }
    } else if (event.source === detachedWindow) {
      if (data.kind === 'mapp-feature-ready') syncDetached();
      else if (data.kind === 'mapp-feature-attach') reattachFeature();
      else if (data.kind === 'mapp-feature-sync' && data.hash === current?.run.hash &&
               typeof data.feature === 'string' && data.feature.length <= 256) {
        selectedFeature = data.feature;
        featureCandidate = 0;
        if (activeTab === 'features') render();
      }
    }
  });

  async function renderOverview(ctx, target, token) {
    const [raw, retained, annotations, differential, blank, qc] = await Promise.all([
      rows(ctx,'matrix_raw',true), rows(ctx,'variable_metadata',true),
      rows(ctx,'annotation_candidates',true), rows(ctx,'differential',true),
      rows(ctx,'blank_filter',true), rows(ctx,'qc_rsd',true)]);
    if (!viewValid(ctx,'overview',token)) return;
    const metadataHeaders = Object.keys(metadata[0]);
    const rawIds = raw?.length ? Object.keys(raw[0]).filter((key) => key !== 'sample') : [];
    const rawCount = raw?.length ? rawIds.length : null;
    const retainedCount = retained?.length ?? null;
    const retainedIds = new Set((retained ?? []).map((row) => row.feature_id));
    const candidates = annotations?.filter((row) => row.source) ?? [];
    const annotated = new Set(candidates.map((row) => row.feature_id));
    const contrasts = [...new Set((differential ?? []).map((row) => row.contrast).filter(Boolean))];
    const retention = [note(rawCount === null || retainedCount === null ?
      'Raw or retained feature table is unavailable in this run.' :
      retainedCount + ' / ' + rawCount + ' features retained (' +
      (100 * retainedCount / Math.max(rawCount,1)).toFixed(1) + '%). ' +
      [...retainedIds].filter((id) => annotated.has(id)).length +
      ' retained features have candidate records. Annotation candidates count records, not features, and can describe features removed by filtering.')];
    if (blank && qc && rawCount !== null && retainedCount !== null) {
      const blankById = new Map(blank.map((row) => [row.feature_id,row]));
      const qcById = new Map(qc.map((row) => [row.feature_id,row]));
      const afterBlank = rawIds.filter((id) => blankById.get(id)?.removed_as_blank === 'FALSE');
      const kept = afterBlank.filter((id) => qcById.get(id)?.removed_by_qc_rsd === 'FALSE');
      const complete = blank.length === rawCount && blankById.size === rawCount &&
        rawIds.every((id) => ['TRUE','FALSE'].includes(blankById.get(id)?.removed_as_blank)) &&
        qc.length === afterBlank.length && qcById.size === afterBlank.length &&
        afterBlank.every((id) => ['TRUE','FALSE'].includes(qcById.get(id)?.removed_by_qc_rsd)) &&
        kept.length === retainedCount && kept.every((id) => retainedIds.has(id));
      if (complete) {
        const settings = ctx.run.preprocessing_filters ?? {};
        const blankSettings = settings.blank_filter ?? {};
        const qcSettings = settings.qc_rsd_filter ?? {};
        retention.push(metrics([
          ['Removed by blank filter',rawCount-afterBlank.length],
          ['Remaining after blank filter',afterBlank.length],
          ['Removed by QC-RSD filter',afterBlank.length-kept.length],
          ['Remaining for analysis',kept.length]
        ]));
        if (Number.isFinite(number(blankSettings.minimum_blank_prevalence)) &&
            Number.isFinite(number(blankSettings.minimum_sample_to_blank_ratio))) {
          retention.push(note('Blank filter: measured in at least ' +
            (100*number(blankSettings.minimum_blank_prevalence)) +
            '% of blank injections, with biological-sample median below ' +
            number(blankSettings.minimum_sample_to_blank_ratio) +
            '× the blank median. Both conditions must hold.'));
        }
        if (Number.isFinite(number(qcSettings.threshold_percent))) {
          retention.push(note('QC-RSD filter: among features passing the blank filter, remove those with QC relative standard deviation above ' +
            number(qcSettings.threshold_percent) + '% (or a non-finite RSD). Calculated across the QC injections after ' +
            (settings.missing_values?.method === 'half_minimum' ? 'half-minimum' :
              label(settings.missing_values?.method)) + ' imputation and ' +
            (settings.normalization?.method === 'pqn' ? 'PQN' :
              label(settings.normalization?.method)) + ' normalization.'));
        }
        retention.push(note('The table below lists each measured feature and its first exclusion stage. Blank-removed features were not tested for QC-RSD; click an ID to inspect its saved raw intensities. Filtering is separate from annotation coverage and downstream PCA/PLS-DA.'));
        const headers = ['Feature ID','Outcome','Sample median','Blank median',
          'Blank prevalence (%)','Sample / blank ratio','QC RSD (%)'];
        const auditRows = rawIds.map((id) => {
          const b = blankById.get(id), q = qcById.get(id);
          return [id,b.removed_as_blank === 'TRUE' ? 'Removed: blank filter' :
            q.removed_by_qc_rsd === 'TRUE' ? 'Removed: QC-RSD' : 'Retained',
          b.sample_median,b.blank_median,(100*number(b.blank_prevalence)).toString(),
          b.sample_to_blank_ratio,q?.qc_rsd_percent ?? ''];
        });
        retention.push(el('details',{class:'retention-audit'},
          el('summary',{text:'Inspect all '+rawCount+' feature filter decisions'}),
          table(headers,auditRows,(id) => openFeature(id),'feature-filter-audit.tsv')));
      } else retention.push(note('Saved blank/QC filter diagnostics do not reconcile with the raw and retained feature IDs; the feature-by-feature audit cannot be shown safely.',true));
    } else retention.push(note('This run has no complete saved blank/QC filter diagnostics; individual exclusion reasons are unavailable.'));
    target.replaceChildren(metrics([
      ['Samples in dataset',metadata.length], ['Raw features',rawCount],
      ['Retained features',retainedCount], ['Annotation candidates',candidates.length], ['Saved contrasts',contrasts.length]
    ]), card('Feature retention & annotation',...retention),
    card('Samples by ' + defaultGroup,
      table(['Group','Samples'], [...new Set(metadata.map((row) => groupName(row)))].sort().map((name) =>
        [name,metadata.filter((row) => groupName(row) === name).length]),null,'sample-groups.tsv')),
    card('Sample roles',table(['Role','Samples'], [...new Set(metadata.map((row) => label(row.sample_type)))].sort().map((role) => [role,metadata.filter((row) => label(row.sample_type) === role).length]),null,'sample-roles.tsv')),
    card('Sample metadata',table(metadataHeaders,metadata.map((row) =>
      metadataHeaders.map((key) => row[key])),null,'sample-metadata.tsv')),
    card('Saved planned contrasts', contrasts.length ? note(contrasts.join(' · ')) : note('No saved differential contrasts in this run.')));
  }

  function plotNode(target, large = false) { const node = el('div', {class:'plot' + (large ? ' large' : ''), role:'img', 'aria-label':'Interactive saved-data plot'}); target.append(node); return node; }
  function plotting(ctx, node, traces, layout) {
    if (!valid(ctx)) return;
    if (!window.Plotly) throw new Error('Local plotly.min.js did not load. Interactive plots are unavailable.');
    Plotly.newPlot(node, traces, {...layout, autosize:true, paper_bgcolor:'transparent', plot_bgcolor:'transparent',
      font:{family:'system-ui, sans-serif',color:'#36505a'}, margin:{l:65,r:20,t:32,b:65,...layout.margin},
      legend:{itemclick:'toggle',itemdoubleclick:'toggleothers',
        ...(traces.some((trace) => trace.marker?.showscale) ? {x:1,y:0,yanchor:'top'} : {})}},
      {responsive:true, displaylogo:false});
  }
  function numericAxes(data, prefix) {
    return data.length ? Object.keys(data[0]).filter((key) => new RegExp('^' + prefix + '[1-9][0-9]*$').test(key)).sort((a,b) => Number(a.slice(prefix.length)) - Number(b.slice(prefix.length))) : [];
  }
  async function renderOrdination(ctx, target, token, stem, prefix) {
    const [scores,colorData] = await Promise.all([
      rows(ctx,stem + '_scores',true),
      stem === 'pca' && available(ctx,'pca_loadings') ? featureColorData(ctx) : null
    ]);
    if (!viewValid(ctx,stem,token)) return;
    if (!scores?.length) { message(target,'No saved ' + stem.toUpperCase() + ' scores for this run.'); return; }
    const axes = numericAxes(scores,prefix);
    if (!axes.length) { message(target,'No saved components found in ' + stem + '_scores.tsv.',true); return; }
    const ui = el('div', {class:'controls'});
    const choices = axes.map((axis) => [axis,axis]);
    if (axes.length >= 3) ui.append(control('View',stem+'-axis-view',[['2D','2D'],['3D','3D']]));
    ui.append(control('X axis',stem+'-axis-x',choices,axes[0]));
    if (axes.length >= 2) ui.append(control('Y axis',stem+'-axis-y',choices,axes[1]));
    if (axes.length >= 3) ui.append(control('Z axis',stem+'-axis-z',choices,axes[2]));
    const wrapper = card('Sample scores',ui,note('Saved scores only. Biological scores may not include QC / blank samples.'));
    const plot = plotNode(wrapper,true);
    const loading = stem === 'pca' && available(ctx,'pca_loadings') ? card('Feature loadings',note('These points represent features, not samples. Click a point to inspect the feature.')) : null;
    const loadingPlot = loading ? plotNode(loading) : null;
    const loadingColors = loading ? featureDisplayOptions(ctx,colorData,'pca',colorData.featureIds,
      loadingPlot,() => redraw().catch((error) => {
        if (viewValid(ctx,stem,token)) message(loadingPlot,error.message,true);
      })) : null;
    if (loadingColors) loading.insertBefore(loadingColors.node,loadingPlot);
    target.replaceChildren(wrapper,...(loading ? [loading] : []),card('Saved original plots',downloadButtons(ctx,savedPlots(ctx,stem))));
    let syncingCamera = false;
    function connectCameras() {
      const relay = (source,destination) => {
        source.cameraRelay = (changes) => {
          const camera = changes['scene.camera'];
          if (!camera || syncingCamera || !viewValid(ctx,stem,token)) return;
          syncingCamera = true;
          Promise.resolve(Plotly.relayout(destination,{'scene.camera':camera}))
            .catch((error) => { if (viewValid(ctx,stem,token)) message(destination,'Could not synchronize 3D rotation: '+error.message,true); })
            .finally(() => { syncingCamera = false; });
        };
        source.on('plotly_relayout',source.cameraRelay);
      };
      relay(plot,loadingPlot);
      relay(loadingPlot,plot);
    }
    async function redraw() {
      for (const node of [plot,loadingPlot]) {
        if (node?.cameraRelay) {
          node.removeListener('plotly_relayout',node.cameraRelay);
          node.cameraRelay = null;
        }
      }
      const axis = (name) => ui.querySelector('#'+stem+'-axis-'+name);
      const view = axis('view')?.value === '3D';
      const x = axis('x').value, y = axis('y')?.value, z = view ? axis('z')?.value : null;
      if (new Set([x,y,z].filter(Boolean)).size !== [x,y,z].filter(Boolean).length) {
        message(plot,'Choose distinct saved components for each axis.',true);
        if (loadingPlot) message(loadingPlot,'Choose distinct saved components for each axis.',true);
        return;
      }
      const selected = new Map(metadata.map((row) => [sampleKey(row),row]));
      const included = scores.filter((row) => selected.has(row.sample) && [x,y,z].filter(Boolean).every((axis) => Number.isFinite(number(row[axis]))));
      if (!included.length) { message(plot,'No finite saved scores for the chosen components.'); return; }
      const names = [...new Set(included.map((row) => groupName(selected.get(row.sample))))].sort();
      const traces = names.map((name,i) => {
        const data = included.filter((row) => groupName(selected.get(row.sample)) === name);
        return {type:view?'scatter3d':'scatter',mode:'markers',name:escapePlot(name),x:data.map((row) => number(row[x])),
          y:data.map((row,j) => y ? number(row[y]) : ((j % 9)-4)*.02),
          ...(view ? {z:data.map((row) => number(row[z]))} : {}),
          customdata:data.map((row) => escapePlot(selected.get(row.sample).sample_id || '(sample identifier unavailable)')),
          hovertemplate:'Sample: %{customdata}<br>' + x + ': %{x:.3f}' + (y ? '<br>' + y + ': %{y:.3f}' : '') + (view ? '<br>' + z + ': %{z:.3f}' : '') + '<extra></extra>',
          marker:{size:view?5:9,color:groupColor(name),opacity:.82},legendgroup:name};
      });
      plot.replaceChildren();
      plot.dataset.sampleGroups = 'true';
      plotting(ctx,plot,traces,view ? {scene:{xaxis:{title:x},yaxis:{title:y},zaxis:{title:z}}} : {xaxis:{title:x},yaxis:{title:y || '1D jitter (display only)'},hovermode:'closest'});
      if (!loadingPlot) return;
      try {
        const [loadings,hoverIndex] = await Promise.all([rows(ctx,'pca_loadings'),featureHoverIndex(ctx)]);
        if (!viewValid(ctx,stem,token) || axis('x')?.value !== x ||
            axis('y')?.value !== y || (axis('view')?.value === '3D') !== view ||
            (view && axis('z')?.value !== z)) return;
        const matching = loadings.filter((row) => row.feature_id && [x,y,z].filter(Boolean).every((axis) => Number.isFinite(number(row[axis]))));
        if (!matching.length) { message(loadingPlot,'No saved loadings for the selected components.'); return; }
        const groups = loadingColors.groups(matching,() =>
          [['Annotated (saved GNPS, ISDB or SIRIUS SMILES)',matching.filter((row) =>
            hoverIndex.features.get(row.feature_id)?.structures.size),'#257e8b',''],
          ['Not annotated (no saved structural SMILES)',matching.filter((row) =>
            !hoverIndex.features.get(row.feature_id)?.structures.size),'#a4a9af','']]);
        const traces = groups.filter(([,data]) => data.length)
          .map(([name,data,color,legendgroup,scale,showlegend]) => ({
          legendgroup,showlegend:showlegend !== false,
          type:view?'scatter3d':'scatter',mode:'markers',name:escapePlot(name),sourceRows:data,
          x:data.map((row) => number(row[x])),
          y:data.map((row,i) => y ? number(row[y]) : (i % 9 - 4)*.02),
          ...(view ? {z:data.map((row) => number(row[z]))} : {}),
          customdata:data.map((row) => escapePlot(row.feature_id)),
          hovertemplate:'Feature: %{customdata}<extra></extra>',
          marker:{size:view?4:7,color,opacity:.7,...scale}
        }));
        loadingPlot.replaceChildren();
        plotting(ctx,loadingPlot,traces.map(({sourceRows,...trace}) => trace),view ?
          {scene:{xaxis:{title:x+' loading'},yaxis:{title:y+' loading'},zaxis:{title:z+' loading'}}} :
          {xaxis:{title:x+' loading'},yaxis:{title:y ? y+' loading' : '1D jitter (display only)'}});
        const featureId = (point) => traces[point.curveNumber]?.sourceRows[point.pointNumber]?.feature_id;
        loadingPlot.on('plotly_click',(event) => {
          const point = event.points[0];
          if (point && valid(ctx)) openFeature(featureId(point));
        });
        bindFeatureHover(ctx,loadingPlot,hoverIndex,featureId);
        if (view) connectCameras();
      } catch (error) { if (valid(ctx)) message(loadingPlot,error.message,true); }
    }
    ui.addEventListener('change',() => redraw().catch((error) => {
      if (viewValid(ctx,stem,token)) message(plot,error.message,true);
    }));
    await redraw();
  }

  async function renderPlsda(ctx, target, token) {
    const statusRows = await rows(ctx,'plsda_status',true);
    if (!viewValid(ctx,'plsda',token)) return;
    if (!statusRows?.length) {
      message(target,'No PLS-DA analysis was saved for this run. Select a completed supervised run.');
      return;
    }
    if (statusRows[0].status !== 'validated') {
      message(target,'PLS-DA '+label(statusRows[0].status)+': '+label(statusRows[0].reason));
      return;
    }
    const [summaries,predictions,permutations,scores,loadings,vip] = await Promise.all([
      rows(ctx,'plsda_summary'),rows(ctx,'plsda_fold_predictions'),
      rows(ctx,'plsda_permutations'),rows(ctx,'plsda_scores'),
      rows(ctx,'plsda_loadings',true),rows(ctx,'plsda_vip',true)
    ]);
    if (!viewValid(ctx,'plsda',token)) return;
    const summary = summaries[0];
    const axes = numericAxes(scores,'LV');
    if (!summary || !predictions.length || !permutations.length || !scores.length || !axes.length ||
        !['balanced_accuracy','accuracy','permutation_p_value','null_mean_balanced_accuracy']
          .every((key) => Number.isFinite(number(summary[key])))) {
      throw new Error('Saved PLS-DA validation outputs are incomplete.');
    }
    const validation = card('Held-out validation',
      metrics([['Balanced accuracy',summary.balanced_accuracy],['Accuracy',summary.accuracy],
        ['Permutation p-value',summary.permutation_p_value],
        ['Null mean balanced accuracy',number(summary.null_mean_balanced_accuracy).toFixed(3)],
        ['Biological samples',summary.samples],['Groups',summary.groups],
        ['Cross-validation',summary.folds+' folds × '+summary.repeats+' repeats'],
        ['Label permutations',summary.permutations]]),
      note('These are saved repeated held-out results for a small dataset, not an independent confirmation. A validated status means the protocol finished, not that discrimination is established.'));
    const ui = el('div',{class:'controls'});
    const choices = axes.map((axis) => [axis,axis]);
    if (axes.length >= 3) ui.append(control('View','plsda-axis-view',[['2D','2D'],['3D','3D']]));
    ui.append(control('X axis','plsda-axis-x',choices,axes[0]));
    if (axes.length >= 2) ui.append(control('Y axis','plsda-axis-y',choices,axes[1]));
    if (axes.length >= 3) ui.append(control('Z axis','plsda-axis-z',choices,axes[2]));
    const scoreCard = card('Full-fit sample scores',ui,
      note('Illustration of the fitted model only; these points are not held-out predictions.'));
    const scorePlot = plotNode(scoreCard,true);
    const loadingCard = card('Full-fit feature loadings',
      note('Illustrative feature contributions to the full-fit model only, not held-out predictive validation. Click a point to inspect its saved feature details.'));
    const loadingPlot = loadings?.length ? plotNode(loadingCard,true) : null;
    let loadingColors, hoverIndex;
    if (loadingPlot) {
      const [colorData,index] = await Promise.all([featureColorData(ctx),featureHoverIndex(ctx)]);
      if (!viewValid(ctx,'plsda',token)) return;
      hoverIndex = index;
      loadingColors = featureDisplayOptions(ctx,colorData,'plsda',loadings.map((row) => row.feature_id).filter(Boolean),
        loadingPlot,redraw);
      loadingCard.insertBefore(loadingColors.node,loadingPlot);
    } else loadingCard.append(note('Feature loadings are unavailable for this earlier run; its saved scores and validation remain available.'));
    const vipCard = card('Full-fit variable importance (VIP)',
      note('Saved full-fit model VIP values describe fitted features; they are not held-out feature validation.'));
    const permutationCard = card('Saved label-permutation distribution',
      note('The vertical line is the saved held-out balanced accuracy. The p-value above comes from the saved validation summary.'));
    const permutationPlot = plotNode(permutationCard);
    let vipPlot, top;
    if (vip?.length) {
      const ranked = vip.filter((row) => row.feature_id && Number.isFinite(number(row.vip)))
        .sort((a,b) => number(b.vip)-number(a.vip));
      if (ranked.length) {
        vipCard.append(note('Top '+Math.min(50,ranked.length)+' saved VIP values, ranked by importance. Click a bar or a feature ID below to inspect its saved details.'));
        vipPlot = plotNode(vipCard,true);
        top = ranked.slice(0,50);
        vipPlot.style.height = Math.max(420,top.length*23+120)+'px';
      } else vipCard.append(note('No finite saved VIP values are available for this run.'));
      vipCard.append(card('All saved VIP values',table(['Feature ID','VIP'],
        vip.map((row) => [row.feature_id,row.vip]),(id) => openFeature(id),'plsda-vip.tsv')));
    } else vipCard.append(note('Feature VIP values are unavailable for this earlier run; its saved scores and validation remain available.'));
    target.replaceChildren(validation,scoreCard,loadingCard,vipCard,permutationCard,
      card('Held-out fold predictions',table(Object.keys(predictions[0]),
        predictions.map((row) => Object.values(row)),null,'plsda-fold-predictions.tsv')),
      card('Saved validation details',table(['Metric','Saved value'],
        Object.entries(summary),null,'plsda-summary.tsv')),
      card('Saved label permutations',table(Object.keys(permutations[0]),
        permutations.map((row) => Object.values(row)),null,'plsda-permutations.tsv')));
    if (vipPlot) {
      plotting(ctx,vipPlot,[{type:'bar',orientation:'h',x:top.map((row) => number(row.vip)),
        y:top.map((_,i) => i),customdata:top.map((row) => escapePlot(row.feature_id)),
        marker:{color:'#197b8a'},
        hovertemplate:'Feature: %{customdata}<br>Saved VIP: %{x:.3f}<extra></extra>'}],
      {xaxis:{title:'Saved full-fit VIP'},yaxis:{title:'Feature ID',autorange:'reversed',
        tickmode:'array',tickvals:top.map((_,i) => i),ticktext:top.map((row) => row.feature_id)},
      margin:{l:140}});
      vipPlot.on('plotly_click',(event) => {
        const point = event.points[0];
        if (point && valid(ctx)) openFeature(top[point.pointNumber]?.feature_id);
      });
    }
    const nullValues = permutations.map((row) => number(row.balanced_accuracy))
      .filter((value) => Number.isFinite(value));
    if (!nullValues.length) message(permutationPlot,'No finite saved permutation accuracies.');
    else plotting(ctx,permutationPlot,[{type:'histogram',x:nullValues,name:'Shuffled labels',
      marker:{color:'#7999a0'},hovertemplate:'Balanced accuracy: %{x:.3f}<br>Permutations: %{y}<extra></extra>'}],
    {xaxis:{title:'Balanced accuracy with shuffled labels',range:[0,1]},yaxis:{title:'Permutation count'},
      shapes:[{type:'line',x0:number(summary.balanced_accuracy),x1:number(summary.balanced_accuracy),
        y0:0,y1:1,yref:'paper',line:{color:'#b13f46',width:3}}]});
    function redraw() {
      if (!viewValid(ctx,'plsda',token)) return;
      const axis = (name) => ui.querySelector('#plsda-axis-'+name);
      const view = axis('view')?.value === '3D';
      const x = axis('x').value, y = axis('y')?.value, z = view ? axis('z')?.value : null;
      if (new Set([x,y,z].filter(Boolean)).size !== [x,y,z].filter(Boolean).length) {
        message(scorePlot,'Choose distinct saved components for each axis.',true);
        if (loadingPlot) message(loadingPlot,'Choose distinct saved components for each axis.',true);
        return;
      }
      const included = scores.filter((row) => [x,y,z].filter(Boolean)
        .every((component) => Number.isFinite(number(row[component]))));
      if (!included.length) { message(scorePlot,'No finite saved scores for the chosen components.'); return; }
      const names = [...new Set(included.map((row) => label(row[ctx.run.design_group])))].sort();
      const sampleIds = new Map(metadata.map((row) => [sampleKey(row),row.sample_id || row.filename]));
      const traces = names.map((name,i) => {
        const data = included.filter((row) => label(row[ctx.run.design_group]) === name);
        return {type:view?'scatter3d':'scatter',mode:'markers',name:escapePlot(name),
          x:data.map((row) => number(row[x])),
          y:data.map((row,j) => y ? number(row[y]) : ((j % 9)-4)*.02),
          ...(view ? {z:data.map((row) => number(row[z]))} : {}),
          customdata:data.map((row) => escapePlot(sampleIds.get(row.sample) || row.sample)),
          hovertemplate:'Sample: %{customdata}<br>'+x+': %{x:.3f}'+
            (y ? '<br>'+y+': %{y:.3f}' : '')+(view ? '<br>'+z+': %{z:.3f}' : '')+'<extra></extra>',
          marker:{size:view?5:9,color:groupColor(name),opacity:.82},legendgroup:name};
      });
      scorePlot.replaceChildren();
      scorePlot.dataset.sampleGroups = 'true';
      plotting(ctx,scorePlot,traces,view ?
        {scene:{xaxis:{title:x},yaxis:{title:y},zaxis:{title:z}}} :
        {xaxis:{title:x},yaxis:{title:y || '1D jitter (display only)'},hovermode:'closest'});
      if (!loadingPlot) return;
      const matching = loadings.filter((row) => row.feature_id && [x,y,z].filter(Boolean)
        .every((component) => Number.isFinite(number(row[component]))));
      if (!matching.length) {
        message(loadingPlot,'No finite saved feature loadings for the chosen components.');
        return;
      }
      const groups = loadingColors.groups(matching,() =>
        [['Annotated (saved GNPS, ISDB or SIRIUS SMILES)',matching.filter((row) =>
          hoverIndex.features.get(row.feature_id)?.structures.size),'#257e8b',''],
        ['Not annotated (no saved structural SMILES)',matching.filter((row) =>
          !hoverIndex.features.get(row.feature_id)?.structures.size),'#a4a9af','']]);
      const loadingTraces = groups.filter(([,data]) => data.length)
        .map(([name,data,color,legendgroup,scale,showlegend]) => ({
        legendgroup,showlegend:showlegend !== false,
        type:view?'scatter3d':'scatter',mode:'markers',name:escapePlot(name),sourceRows:data,
        x:data.map((row) => number(row[x])),
        y:data.map((row,i) => y ? number(row[y]) : ((i % 9)-4)*.02),
        ...(view ? {z:data.map((row) => number(row[z]))} : {}),
        customdata:data.map((row) => escapePlot(row.feature_id)),
        hovertemplate:'Feature: %{customdata}<br>'+x+': %{x:.3f}'+
          (y ? '<br>'+y+': %{y:.3f}' : '')+(view ? '<br>'+z+': %{z:.3f}' : '')+'<extra></extra>',
        marker:{size:view?4:7,color,opacity:.7,...scale}
      }));
      hideFeatureHover();
      if (loadingPlot.data && window.Plotly) Plotly.purge(loadingPlot);
      loadingPlot.replaceChildren();
      plotting(ctx,loadingPlot,loadingTraces.map(({sourceRows,...trace}) => trace),view ?
        {scene:{xaxis:{title:x+' loading'},yaxis:{title:y+' loading'},zaxis:{title:z+' loading'}}} :
        {xaxis:{title:x+' loading'},yaxis:{title:y ? y+' loading' : '1D jitter (display only)'},hovermode:'closest'});
      const featureId = (point) => loadingTraces[point.curveNumber]?.sourceRows[point.pointNumber]?.feature_id;
      loadingPlot.on('plotly_click',(event) => {
        const point = event.points[0];
        if (point && valid(ctx)) openFeature(featureId(point));
      });
      bindFeatureHover(ctx,loadingPlot,hoverIndex,featureId);
    }
    ui.addEventListener('change',redraw);
    redraw();
  }

  async function renderVolcano(ctx,target,token) {
    const [data,hoverIndex,colorData] = await Promise.all([
      rows(ctx,'differential',true),featureHoverIndex(ctx),featureColorData(ctx)
    ]);
    if (!viewValid(ctx,'volcano',token)) return;
    if (!data?.length) { message(target,'No saved differential tests for this run.'); return; }
    const comparisons = new Map();
    for (const row of data) {
      if (!row.contrast || !row.numerator || !row.denominator || row.numerator === row.denominator) continue;
      if (!comparisons.has(row.contrast))
        comparisons.set(row.contrast,{name:row.contrast,numerator:row.numerator,denominator:row.denominator,rows:[]});
      comparisons.get(row.contrast).rows.push(row);
    }
    if (!comparisons.size) { message(target,'No named group comparisons found in differential.tsv.',true); return; }
    const all = [...comparisons.values()];
    const groups = [...new Set([...metadata.map((row) => row[ctx.run.design_group]),
      ...all.flatMap((item) => [item.numerator,item.denominator])].filter((name) => name && name !== 'NA'))].sort();
    const validGroups = new Set(groups);
    let left = all[0].denominator, right = all[0].numerator, selectedGroup = '';
    const pool = el('div',{class:'volcano-groups'});
    const instruction = note('Drag a group card to Left or Right, or select a card then activate a target. The horizontal effect is Right relative to Left.');
    const chooser = el('div',{class:'volcano-chooser'},pool);
    const sides = {};
    function assign(side,group) {
      if (!validGroups.has(group)) return;
      if (side === 'left') {
        if (group === right) right = left;
        left = group;
      } else {
        if (group === left) left = right;
        right = group;
      }
      update();
    }
    for (const [side,title] of [['left','Left · reference'],['right','Right · numerator']]) {
      const slot = el('button',{type:'button',class:'volcano-target',draggable:'true'});
      slot.addEventListener('click',() => { if (selectedGroup) assign(side,selectedGroup); });
      slot.addEventListener('dragstart',(event) => {
        event.dataTransfer.setData('text/plain',side === 'left' ? left : right);
        event.dataTransfer.effectAllowed = 'copyMove';
      });
      slot.addEventListener('dragover',(event) => { event.preventDefault(); event.dataTransfer.dropEffect = 'move'; slot.classList.add('drop-ready'); });
      slot.addEventListener('dragleave',() => slot.classList.remove('drop-ready'));
      slot.addEventListener('drop',(event) => {
        event.preventDefault();
        slot.classList.remove('drop-ready');
        assign(side,event.dataTransfer.getData('text/plain'));
      });
      sides[side] = slot;
      chooser.append(el('div',{class:'volcano-side'},el('strong',{text:title}),slot));
    }
    const swap = el('button',{type:'button',class:'button secondary small',text:'Swap sides',onClick:() => {
      [left,right] = [right,left]; update();
    }});
    chooser.append(swap);
    const description = note('');
    const cardPlot = card('Saved differential tests',instruction,chooser,description,
      note('q ≤ 0.05 marks saved adjusted significance. Filters on samples do not alter p or q.'));
    const plot = plotNode(cardPlot,true);
    const colorOptions = featureDisplayOptions(ctx,colorData,'volcano',
      [...new Set(data.map((row) => row.feature_id).filter(Boolean))],plot,update);
    const significanceHint = note('With annotation colours, diamonds mark saved q ≤ 0.05; circles mark other or missing q-values.');
    significanceHint.hidden = true;
    cardPlot.insertBefore(colorOptions.node,plot);
    cardPlot.insertBefore(significanceHint,plot);
    const links = card('Saved original plots');
    target.replaceChildren(cardPlot,links);
    const cards = groups.map((group) => {
      const button = el('button',{type:'button',class:'volcano-group',draggable:'true',text:group,'aria-pressed':'false'});
      button.addEventListener('click',() => { selectedGroup = group; update(); });
      button.addEventListener('dragstart',(event) => {
        selectedGroup = group;
        event.dataTransfer.setData('text/plain',group);
        event.dataTransfer.effectAllowed = 'copyMove';
        update();
      });
      pool.append(button);
      return button;
    });
    function update() {
      for (const [side,group] of [['left',left],['right',right]]) {
        sides[side].textContent = group;
        sides[side].setAttribute('aria-label',(side === 'left' ? 'Left reference' : 'Right numerator') +
          ': '+group+'. Drop a group here or select a card and activate this target.');
      }
      cards.forEach((button,i) => {
        button.classList.toggle('selected',groups[i] === selectedGroup);
        button.classList.toggle('in-use',groups[i] === left || groups[i] === right);
        button.setAttribute('aria-pressed',String(groups[i] === selectedGroup));
      });
      const match = all.find((item) => item.numerator === right && item.denominator === left) ||
        all.find((item) => item.numerator === left && item.denominator === right);
      hideFeatureHover();
      if (plot.data && window.Plotly) Plotly.purge(plot);
      plot.replaceChildren();
      if (!match) {
        description.textContent = right+' vs '+left+' has no saved differential comparison.';
        message(plot,'No saved differential test for this group pair. Choose a pair with a saved comparison.');
        links.replaceChildren(el('h3',{text:'Saved original plots'}),note('No original plot for this group pair.'));
        return;
      }
      const reversed = match.numerator === left;
      description.textContent = 'Viewing '+right+' vs '+left+' from saved contrast '+match.name+
        ' ('+match.numerator+' vs '+match.denominator+'). '+
        (reversed ? 'Saved effect signs are reversed for this display; p and q remain saved values.' :
          'Saved effect signs, p and q are unchanged.');
      significanceHint.hidden = !colorOptions.field();
      const plotted = match.rows.filter((row) => Number.isFinite(number(row.effect)) &&
        number(row.p_value) >= 0 && number(row.p_value) <= 1);
      const finitePositive = plotted.map((row) => number(row.p_value)).filter((value) => value > 0);
      const floor = Math.min(1e-300,...finitePositive);
      const colorGroups = colorOptions.groups(plotted,() =>
        [[false,plotted.filter((row) => !(number(row.q_value) <= .05)),'#83a3a6',''],
          [true,plotted.filter((row) => number(row.q_value) <= .05),'#c35e44','']]
          .map(([significant,subset,color]) =>
            [significant?'q ≤ 0.05':'Other / q missing',subset,color,'']));
      const traces = colorGroups.map(([name,subset,color,legendgroup,scale,showlegend]) => ({
        type:'scattergl',mode:'markers',name:escapePlot(name),legendgroup,
        showlegend:showlegend !== false,
        x:subset.map((row) => (reversed ? -1 : 1)*number(row.effect)),
        y:subset.map((row) => -Math.log10(Math.max(number(row.p_value),floor))),
        customdata:subset.map((row) => [escapePlot(row.feature_id),escapePlot(row.q_value)]),
        hovertemplate:'Feature: %{customdata[0]}<br>Effect: %{x:.3f}<br>−log10(p): %{y:.3f}<br>Saved q: %{customdata[1]}<extra></extra>',
        marker:{color,size:7,opacity:.75,...scale,...(colorOptions.field() ?
          {symbol:subset.map((row) => number(row.q_value) <= .05 ? 'diamond' : 'circle')} : {})},
        sourceRows:subset
      }));
      if (!plotted.length) message(plot,'No finite saved effects and p-values for this comparison.');
      else {
        const extent = plotted.reduce((max,row) => Math.max(max,Math.abs(number(row.effect))),0) || 1;
        const bound = extent * 1.05;
        plotting(ctx,plot,traces.map(({sourceRows,...trace}) => trace),{
          xaxis:{title:'Saved '+escapePlot(label(plotted[0].effect_scale))+' ('+escapePlot(right)+' vs '+escapePlot(left)+')',
            range:[-bound,bound],zeroline:true},
          yaxis:{title:'−log10(saved p-value)'},hovermode:'closest'});
        plot.on('plotly_click',(event) => {
          const p = event.points[0];
          if (p && valid(ctx)) openFeature(traces[p.curveNumber]?.sourceRows[p.pointNumber]?.feature_id);
        });
        bindFeatureHover(ctx,plot,hoverIndex,(point) =>
          traces[point.curveNumber]?.sourceRows[point.pointNumber]?.feature_id);
      }
      links.replaceChildren(el('h3',{text:'Saved original plots'}),
        note('Original '+match.numerator+' vs '+match.denominator+' plots retain their saved orientation.'),
        downloadButtons(ctx,savedPlots(ctx,'volcano_'+match.name)));
    }
    update();
  }

  function openFeature(id) {
    if (!id || !current) return;
    selectedFeature = text(id);
    featureCandidate = 0;
    if (detachedMode) {
      publishFeatureSelection();
      render();
      return;
    }
    if (syncDetached()) {
      closeFeatureDrawer();
      if (activeTab === 'features') render();
      return;
    }
    const drawer = $('feature-drawer');
    const content = $('feature-drawer-content');
    const ctx = current, token = ctx.viewToken, version = ++drawerVersion;
    $('feature-drawer-title').textContent = 'Feature '+selectedFeature+' · details';
    drawer.hidden = false;
    drawer.open = true;
    $('dashboard-main').classList.add('drawer-open');
    message(content,'Loading verified saved feature details…');
    resizeActivePlots();
    renderFeatures(ctx,content,token,true,version).catch((error) => {
      if (viewValid(ctx,activeTab,token) && version === drawerVersion) message(content,error.message || text(error),true);
    });
  }
  function selectedFilterRows(features) {
    const source = $('candidate-source').value;
    const agreement = $('candidate-consensus')?.value || '';
    const npc = ['npc_pathway','npc_superclass','npc_class'].map((field) => $(field).value);
    const taxon = $('taxon-filter').value.trim().toLocaleLowerCase();
    const query = $('candidate-search').value.trim().toLocaleLowerCase();
    return features.filter(({id, records, agreements}) =>
      (!agreement || agreements.has(agreement)) &&
      (!source || records.some((row) => row.source === source)) &&
      npc.every((value,i) => !value || records.some((row) => row[['npc_pathway','npc_superclass','npc_class'][i]] === value)) &&
      (!taxon || records.some((row) => row.source === 'met_annot_enhancer' &&
        text(row.organism).toLocaleLowerCase().includes(taxon))) &&
      (!query || id.toLocaleLowerCase().includes(query) ||
        records.some((row) => Object.values(row).some((value) => text(value).toLocaleLowerCase().includes(query)))));
  }
  async function renderFeatures(ctx,target,token,drawer = false,version = 0) {
    const [candidates, consensus, variables, raw, hoverIndex] = await Promise.all([
      rows(ctx,'annotation_candidates',true),rows(ctx,'annotation_horizontal',true),
      rows(ctx,'variable_metadata',true),rows(ctx,'matrix_raw',true),featureHoverIndex(ctx)]);
    const hostTab = drawer ? activeTab : 'features';
    if (!viewValid(ctx,hostTab,token) ||
        (drawer && !detachedMode && (version !== drawerVersion || !$('feature-drawer').open))) return;
    const ids = new Set([...(candidates ?? []).map((row) => row.feature_id),
      ...(consensus ?? []).map((row) => row.feature_id),
      ...(variables ?? []).map((row) => row.feature_id),
      ...Object.keys(raw?.[0] ?? {}).filter((key) => key !== 'sample')].filter(Boolean));
    const focus = el('div');
    let features, picker, chooser, browse, filterControls, taxonInput, search, results;
    if (!drawer) {
      const byFeature = new Map();
      function feature(id) {
        if (!byFeature.has(id)) byFeature.set(id,{id,records:[],agreements:new Set(),agreeingSources:new Set()});
        return byFeature.get(id);
      }
      for (const row of candidates ?? []) if (row.feature_id) feature(row.feature_id).records.push(row);
      for (const row of consensus ?? []) if (row.feature_id) {
        const item = feature(row.feature_id);
        if (row.sources_number_IK2D !== '') item.agreements.add(row.sources_number_IK2D);
        if (row.sources_IK2D) item.agreeingSources.add(row.sources_IK2D);
      }
      features = [...byFeature.values()].sort((a,b) => a.id.localeCompare(b.id,undefined,{numeric:true}));
      picker = el('input', {type:'search', id:'feature-id', list:'feature-ids', placeholder:'e.g. 467', value:selectedFeature});
      const options = el('datalist',{id:'feature-ids'},...[...ids].map((id) => el('option',{value:id})));
      chooser = card('Inspect a feature', el('label',{for:'feature-id',text:'Feature ID'}),picker,options,
        note('Enter a saved feature ID or select one from the feature table / plots.'));
      browse = card('Search saved features');
      filterControls = el('div',{class:'controls'});
      const sources = [...new Set((candidates ?? []).map((row) => row.source).filter(Boolean))].sort();
      filterControls.append(control('Candidate record source (not structure agreement)','candidate-source',
        [['','All record sources'],...sources.map((v) => [v,v])]));
      if (consensus) filterControls.append(control('Saved InChIKey2D agreement','candidate-consensus',
        [['','All'],...['0','1','2','3'].map((v) => [v,v+' agreeing source'+(v === '1' ? '' : 's')])]));
      else filterControls.append(note('No saved InChIKey2D consensus table in this run.'));
      for (const field of ['npc_pathway','npc_superclass','npc_class']) {
        const values = [...new Set((candidates ?? []).map((row) => row[field]).filter(Boolean))].sort();
        filterControls.append(control(field.replace('_',' ').toUpperCase(),field,[['','All'],...values.map((v) => [v,v])]));
      }
      taxonInput = el('input',{type:'search',id:'taxon-filter',placeholder:'Case-insensitive taxon name'});
      filterControls.append(el('div',{},el('label',{for:'taxon-filter',text:'Reported taxon contains'}),taxonInput));
      search = el('input',{type:'search',id:'candidate-search',placeholder:'ID, source key, label, NPC, taxon…'});
      results = el('div');
      browse.append(filterControls,
        note('Record source matches an annotation candidate row: SIRIUS, CANOPUS class predictions, or met_annot_enhancer (ISDB/MS1 evidence). GNPS can contribute to saved 2D-key agreement without a candidate row here. Selecting a record source does not recalculate or narrow the saved consensus count.'),
        el('label',{for:'candidate-search',text:'Search all saved candidate fields'}),search,
        note('Each feature appears once. Saved InChIKey2D agreement and agreeing names come only from annotation_horizontal.tsv (GNPS, ISDB, SIRIUS). CANOPUS reports classes, not structure consensus; its source-specific record remains available when inspecting a feature.'),results);
      target.replaceChildren(chooser,focus,browse);
    } else target.replaceChildren(focus);
    async function focusFeature() {
      focus.replaceChildren();
      if (!selectedFeature) { focus.append(note('Choose a feature to see its saved intensities, candidate evidence and structures.')); return; }
      if (!ids.has(selectedFeature)) { focus.append(note('Feature ID not found in the saved tables for this run.',true)); return; }
      const id = selectedFeature;
      const mine = (candidates ?? []).filter((row) => row.feature_id === id);
      const counts = (consensus ?? []).filter((row) => row.feature_id === id);
      const stageNames = Object.keys(ctx.run.outputs).filter((path) => /^tables\/matrix_[A-Za-z0-9_.-]+\.tsv$/.test(path))
        .map((path) => path.slice(14,-4));
      const facetFields = sampleGroupColumns.filter((field) => {
        const count = new Set(metadata.map((row) => groupName(row,field))).size;
        return count > 1 && count <= 12;
      });
      const settings = ctx.featurePlotSettings ??= savedFeaturePlotSettings(ctx,stageNames,facetFields);
      const plotCard = card('Descriptive sample intensities');
      const stageSelect = control('Saved intensity stage','feature-stage',stageNames.map((name) => [name,name]),settings.stage);
      const groupControl = control('Group by','intensity-group',sampleGroupColumns.map((name) => [name,name]),settings.group);
      const roleControl = control('Sample role','intensity-role',
        [['','All roles'],...[...new Set(metadata.map((row) => row.sample_type).filter(Boolean))].sort().map((name) => [name,name])],settings.role);
      const facetControl = control('Facet by','intensity-facet',
        [['','No facets'],...facetFields.map((name) => [name,name])],settings.facet);
      const groupsControl = control('Groups to show','intensity-groups',[]);
      const groupSelect = groupControl.querySelector('select');
      const roleSelect = roleControl.querySelector('select');
      const facetSelect = facetControl.querySelector('select');
      const groupsSelect = groupsControl.querySelector('select');
      groupsSelect.multiple = true;
      function resetGroups() {
        const names = [...new Set(metadata.map((row) => groupName(row,groupSelect.value)))].sort();
        selector(groupsSelect,names);
        groupsSelect.size = Math.min(6,Math.max(2,names.length));
        for (const option of groupsSelect.options) option.selected = settings.groups === null || settings.groups.has(option.value);
      }
      resetGroups();
      const displayOptions = el('details',{class:'intensity-options'},el('summary',{text:'Display options'}),
        el('div',{class:'intensity-options-body'},
          el('div',{class:'controls'},stageSelect,groupControl,roleControl,facetControl),
          groupsControl,note('Ctrl/Cmd-click to show multiple groups. These controls affect only this descriptive plot. Group colours are shared in the sidebar.')));
      displayOptions.open = settings.open;
      displayOptions.addEventListener('toggle',() => {
        settings.open = displayOptions.open;
        saveFeaturePlotSettings(ctx);
      });
      plotCard.append(displayOptions,
        note('Saved peak '+label(ctx.run.intensity_measure)+'. Boxes and points are descriptive; no test is recomputed.'));
      const box = plotNode(plotCard);
      const annotation = card('Chemical depictions');
      const structures = galleryStructures(mine,counts);
      const gallery = el('div',{class:'gallery'});
      for (const {ik2d,records} of structures.values()) {
        const names = [...new Set(records.map((row) => text(row.label).trim()).filter((name) =>
          name && name !== 'Saved horizontal structure'))];
        const sources = preferredStructureSources.filter((source) =>
          ik2d && records.some((row) => row.source === source && row.ik2d === ik2d));
        const structure = el('div',{class:'structure'},
          el('strong',{text:names.join(' · ') || 'Saved candidate structure'}));
        if (sources.length) {
          structure.append(el('div',{class:'structure-sources','aria-label':
            'Saved sources sharing 2D key '+ik2d+': '+sources.map((source) => source.toUpperCase()).join(', ')},
          ...sources.map((source) => sourceStar(source,records))));
        } else structure.append(note('No GNPS/ISDB/SIRIUS star: no eligible saved source has a valid 2D key.'));
        if (ik2d) structure.append(el('p',{class:'structure-key',text:'Saved 2D key: '+ik2d}));
        const depiction = records.find((row) => usableSmiles(row.smiles));
        if (depiction) {
          const image = el('img',{class:'structure-depiction',loading:'lazy',referrerpolicy:'no-referrer',
            src:depictionUrl(depiction.smiles,440,300),
            alt:'Representative 2D depiction from saved '+depiction.provenance+' SMILES'+
              (ik2d ? ' for 2D key '+ik2d : '')+'; candidate structure, not confirmed identity'});
          const fallback = el('p',{class:'note',text:'Depiction unavailable; saved source SMILES are below.',hidden:''});
          image.addEventListener('error',() => { image.hidden = true; fallback.hidden = false; });
          structure.append(image,fallback);
        } else structure.append(note('Saved SMILES exceeds the depiction service’s 5000-byte limit.'));
        if (new Set(records.map((row) => row.smiles)).size > 1)
          structure.append(note('One representative SMILES is depicted; saved source SMILES differ.'));
        const details = el('details',{},el('summary',{text:'Saved source SMILES and evidence'}));
        for (const record of records) {
          const score = recordedText(record.score), taxon = recordedText(record.taxon);
          details.append(el('div',{class:'structure-record'},
            el('span',{class:'structure-record-label',text:record.provenance+
              (record.label && record.label !== 'Saved horizontal structure' ? ' · '+record.label : '')+
              (record.ik2d ? ' · 2D key '+record.ik2d : ' · 2D key not recorded')}),
            el('code',{text:record.smiles}),
            score || taxon ? el('span',{class:'structure-record-evidence',text:
              [score ? 'Saved taxonomic score: '+score : '',taxon ? 'Reported taxon: '+taxon : '']
                .filter(Boolean).join(' · ')}) : null));
        }
        structure.append(details);
        gallery.append(structure);
      }
      annotation.append(el('p',{class:'note'},
        document.createTextNode('Candidate structures only. Drawings by the '),
        el('a',{href:'https://api.naturalproducts.net/latest/docs#/depict/depict_2d_molecule_depict_2D_get',
          target:'_blank',rel:'noopener noreferrer',referrerpolicy:'no-referrer',text:'Natural Products API'}),
        document.createTextNode('; each image request sends saved SMILES and your network address.')),
        note('Stars identify saved sources sharing this valid 2D key, not confirmed identity. ISDB star color uses the recorded taxonomic score on a 0–8 visual scale; a neutral star means no score was saved.'),
        structures.size ? gallery : note('No saved SMILES for this feature.'));
      const consensusText = consensus ? (counts.length ? 'Saved InChIKey2D agreement (count and sources): ' +
        counts.map((row) => label(row.sources_number_IK2D)+' ['+label(row.sources_IK2D)+']').join(' · ')
        : 'No saved consensus row for this feature.') : 'No saved consensus table for this run.';
      const evidence = card('Source-specific candidate records for feature '+id, note(consensusText),
        note('Source record keys identify input rows, not other feature IDs. CANOPUS class records are not structure consensus.'));
      if (mine.length) {
        const selected = control('Source-specific evidence','feature-evidence',
          mine.map((row,i) => [text(i),(row.source || 'Source unknown')+' · '+(row.candidate_id || 'key missing')+' · '+(row.label || 'unlabelled')]),
          text(Math.min(featureCandidate,mine.length-1)));
        const details = el('div');
        const updateEvidence = () => {
          featureCandidate = Number(selected.querySelector('select').value);
          const row = mine[featureCandidate];
          details.replaceChildren(table(['Saved field','Value'],Object.entries(row).filter(([,value]) => value !== '').map(([key,value]) => [key,value]),null,'candidate-evidence.tsv'));
        };
        selected.addEventListener('change',updateEvidence);
        evidence.append(selected,details);
        updateEvidence();
      }
      const vm = (variables ?? []).find((row) => row.feature_id === id);
      if (vm) evidence.append(el('h4',{text:'Measured feature metadata'}),
        table(['Saved field','Value'],Object.entries(vm).filter(([key,value]) => key && value).map(([key,value]) => [key,value]),null,'feature-metadata.tsv'));
      const split = el('div',{class:'split'},annotation,plotCard);
      focus.append(split,evidence);
      async function drawBox() {
        const stage = stageSelect.querySelector('select').value;
        const matrix = stage === 'raw' ? raw : await rows(ctx,'matrix_'+stage);
        if (!viewValid(ctx,hostTab,token) ||
            (drawer && !detachedMode && (version !== drawerVersion || !$('feature-drawer').open)) ||
            selectedFeature !== id || stageSelect.querySelector('select').value !== stage) return;
        const fallback = !matrix?.length || !Object.hasOwn(matrix[0],id);
        if (fallback && (!raw?.length || !Object.hasOwn(raw[0],id))) { message(box,'Feature absent from selected stage and raw matrix.',true); return; }
        const actual = fallback ? raw : matrix;
        const groupField = groupSelect.value, facetField = facetSelect.value;
        const allowed = new Set([...groupsSelect.selectedOptions].map((option) => option.value));
        const selected = new Map(metadata.filter((row) =>
          (!roleSelect.value || row.sample_type === roleSelect.value) &&
          allowed.has(groupName(row,groupField))).map((row) => [sampleKey(row),row]));
        const values = actual.filter((row) => selected.has(row.sample) && Number.isFinite(number(row[id])));
        hideFeatureHover();
        if (box.data && window.Plotly) Plotly.purge(box);
        box.replaceChildren();
        split.classList.remove('faceted');
        if (!values.length) {
          message(box,'No finite intensities match the controls in this plot.');
          return;
        }
        const names = [...new Set(values.map((row) => groupName(selected.get(row.sample),groupField)))].sort();
        const facetNames = facetField ?
          [...new Set(values.map((row) => groupName(selected.get(row.sample),facetField)))].sort() : [''];
        const faceted = facetNames.length > 1;
        split.classList.toggle('faceted',faceted);
        const columns = faceted ? Math.min(2,facetNames.length) : 1;
        const rowsCount = Math.ceil(facetNames.length / columns);
        box.style.height = faceted ? (rowsCount*290+90)+'px' : '';
        const layout = {hovermode:'closest'};
        if (faceted) layout.annotations = [];
        const traces = [], shown = new Set();
        for (const [index,facet] of facetNames.entries()) {
          const suffix = index ? String(index+1) : '';
          if (faceted) {
            const col = index % columns, row = Math.floor(index / columns);
            const gapX = .08, gapY = .12;
            const width = (1-gapX*(columns-1))/columns;
            const height = (1-gapY*(rowsCount-1))/rowsCount;
            const left = col*(width+gapX), top = 1-row*(height+gapY);
            layout['xaxis'+suffix] = {domain:[left,left+width],anchor:'y'+suffix,
              categoryorder:'array',categoryarray:names,title:row === rowsCount-1 ? 'Descriptive group' : ''};
            layout['yaxis'+suffix] = {domain:[top-height,top],anchor:'x'+suffix,
              ...(index ? {matches:'y'} : {}),title:col === 0 ? 'Saved '+(fallback?'raw':stage)+' intensity' : ''};
            layout.annotations.push({x:left+width/2,y:Math.min(1.04,top+.025),xref:'paper',yref:'paper',
              text:escapePlot(facetField+': '+facet),showarrow:false,font:{size:12,color:'#244d56'}});
          }
          for (const name of names) {
            const subset = values.filter((row) =>
              groupName(selected.get(row.sample),groupField) === name &&
              (!facetField || groupName(selected.get(row.sample),facetField) === facet));
            if (!subset.length) continue;
            const color = groupColor(name);
            traces.push({type:'box',name:escapePlot(name),legendgroup:name,showlegend:!shown.has(name),
              ...(faceted ? {xaxis:'x'+suffix,yaxis:'y'+suffix} : {}),
              y:subset.map((row) => number(row[id])),
              customdata:subset.map((row) => escapePlot(selected.get(row.sample).sample_id || '(sample identifier unavailable)')),
              boxpoints:'all',jitter:.3,pointpos:0,marker:{color,size:5,opacity:.7},line:{color},
              hovertemplate:'Sample: %{customdata}<br>Intensity: %{y:.3f}<extra></extra>'});
            shown.add(name);
          }
        }
        if (!faceted) {
          layout.xaxis = {title:'Descriptive group'};
          layout.yaxis = {title:'Saved '+(fallback?'raw':stage)+' intensity'};
        }
        box.dataset.sampleGroups = 'true';
        plotting(ctx,box,traces,layout);
        bindFeatureHover(ctx,box,hoverIndex,() => id);
        const old = plotCard.querySelector('.fallback');
        if (old) old.remove();
        if (fallback) { const warning = note('Feature is absent from '+stage+' matrix; showing saved raw intensities instead.'); warning.classList.add('fallback'); plotCard.append(warning); }
      }
      const refresh = () => drawBox().catch((error) => {
        if (viewValid(ctx,hostTab,token) && (!drawer || version === drawerVersion)) message(box,error.message,true);
      });
      stageSelect.addEventListener('change',() => {
        settings.stage = stageSelect.querySelector('select').value;
        saveFeaturePlotSettings(ctx);
        refresh();
      });
      groupControl.addEventListener('change',() => {
        settings.group = groupSelect.value;
        settings.groups = null;
        resetGroups();
        saveFeaturePlotSettings(ctx);
        refresh();
      });
      roleControl.addEventListener('change',() => {
        settings.role = roleSelect.value;
        saveFeaturePlotSettings(ctx);
        refresh();
      });
      facetControl.addEventListener('change',() => {
        settings.facet = facetSelect.value;
        saveFeaturePlotSettings(ctx);
        refresh();
      });
      groupsControl.addEventListener('change',() => {
        settings.groups = new Set([...groupsSelect.selectedOptions].map((option) => option.value));
        saveFeaturePlotSettings(ctx);
        refresh();
      });
      await drawBox();
    }
    function browseRows() {
      const filtered = selectedFilterRows(features);
      const joined = (records,field) => [...new Set(records.map((row) => row[field]).filter((value) => value && value !== 'NA'))].join(' · ');
      const t = table(['Feature ID','Saved InChIKey2D agreement','Saved agreeing structure sources','Source record keys','Candidate labels','NPC classes','Reported taxa'],
        filtered.map(({id,records,agreements,agreeingSources}) => [id,[...agreements].sort().join(' · '),
          [...agreeingSources].join(' · '),joined(records,'candidate_id'),joined(records,'label'),
          joined(records,'npc_class'),joined(records,'organism')]),(id) => {
          selectedFeature = id; picker.value = id; featureCandidate = 0;
          publishFeatureSelection();
          focusFeature().catch((error) => { if (valid(ctx)) message(focus,error.message,true); });
          chooser.scrollIntoView({behavior:'smooth',block:'start'});
        },'feature-browse.tsv');
      results.replaceChildren(note(filtered.length+' features match the feature search and selectors.'),t);
    }
    if (!drawer) {
      picker.addEventListener('change',() => {
        selectedFeature = picker.value.trim(); featureCandidate = 0;
        publishFeatureSelection();
        focusFeature().catch((error) => { if (valid(ctx)) message(focus,error.message,true); });
      });
      search.addEventListener('input',browseRows);
      taxonInput.addEventListener('input',browseRows);
      filterControls.addEventListener('change',browseRows);
      browseRows();
    }
    await focusFeature();
  }

  async function renderClasses(ctx,target,token) {
    const [candidates,raw] = await Promise.all([rows(ctx,'annotation_candidates',true),rows(ctx,'matrix_raw',true)]);
    if (!viewValid(ctx,'classes',token)) return;
    if (!candidates?.length || !raw?.length) {
      message(target,'This sealed run has no saved CANOPUS candidates or raw intensity matrix for chemical-class exploration.');
      return;
    }
    const biological = metadata.filter((row) => row.sample_type === 'sample');
    const rawBySample = new Map(raw.map((row) => [row.sample,row]));
    if (!biological.length || biological.some((row) => !rawBySample.has(sampleKey(row)))) {
      message(target,'Saved biological-sample metadata and raw intensity rows do not match; no comparison can be shown.',true);
      return;
    }
    const groupFields = sampleGroupColumns.filter((field) =>
      new Set(biological.map((row) => row[field]).filter((value) => value && value !== 'NA')).size >= 2);
    if (!groupFields.length) {
      message(target,'No saved categorical biological grouping field has two groups in this run.');
      return;
    }
    const predicted = new Map();
    for (const row of candidates) {
      if (row.source !== 'canopus' || !row.feature_id || !Object.hasOwn(raw[0],row.feature_id)) continue;
      const previous = predicted.get(row.feature_id);
      if (!previous || number(row.class_probability) > number(previous.class_probability)) predicted.set(row.feature_id,row);
    }
    if (!predicted.size) {
      message(target,'No CANOPUS-predicted features have saved raw measurements in this run.');
      return;
    }
    const root = {name:'All NPC pathways',level:0,ids:[],children:new Map()};
    const taxonomy = ['npc_pathway','npc_superclass','npc_class'];
    for (const [id,row] of predicted) {
      let node = root;
      node.ids.push(id);
      for (const field of taxonomy) {
        const name = row[field] && row[field] !== 'NA' ? row[field] : 'Unclassified';
        if (!node.children.has(name)) node.children.set(name,{name,level:node.level+1,ids:[],children:new Map()});
        node = node.children.get(name);
        node.ids.push(id);
      }
      node.children.set(id,{name:id,level:4,ids:[id],children:new Map()});
    }
    const state = ctx.classView ??= {field:groupFields.includes(defaultGroup) ? defaultGroup : groupFields[0],
      first:'',second:'',path:[]};
    if (!groupFields.includes(state.field)) state.field = groupFields[0];
    const field = control('Compare biological groups by','class-group-field',groupFields.map((value) => [value,value]),state.field);
    const first = control('Group A','class-group-a',[]);
    const second = control('Group B','class-group-b',[]);
    const fieldSelect = field.querySelector('select'), aSelect = first.querySelector('select'), bSelect = second.querySelector('select');
    const controls = el('div',{class:'controls'},field,first,second);
    const context = note('CANOPUS class assignments are predictions, not confirmed chemical identities. Each feature contributes once to one saved pathway → superclass → class branch. Raw peak heights are summed per biological sample, then averaged within each selected group; no statistical tests or normalization are performed.');
    const summary = el('p',{class:'notice','aria-live':'polite'});
    const trail = el('nav',{class:'class-trail','aria-label':'NPC hierarchy'});
    const results = el('div');
    const chart = plotNode(card('Group composition · share of parent'));
    const chartCard = chart.parentNode;
    chartCard.append(note('Each bar is the percentage of the current parent’s mean summed raw intensity contributed by that child. The table shows absolute group means and B/A; unequal sample counts are averaged per sample.'));
    target.replaceChildren(card('Compare saved chemical classes',controls,context,summary),card('Explore pathway → superclass → class → feature',trail,results),chartCard);
    function groupOptions() {
      const values = [...new Set(biological.map((row) => row[state.field]).filter((value) =>
        value && value !== 'NA'))].sort((a,b) => a.localeCompare(b));
      if (!values.includes(state.first)) state.first = values[0];
      if (!values.includes(state.second) || state.second === state.first) state.second = values.find((value) => value !== state.first);
      selector(aSelect,values);
      selector(bSelect,values);
      aSelect.value = state.first;
      bSelect.value = state.second;
    }
    function branch() {
      let node = root;
      const path = [];
      for (const name of state.path) {
        if (!node.children.has(name)) break;
        node = node.children.get(name);
        path.push(name);
      }
      state.path = path;
      return node;
    }
    function update() {
      if (!viewValid(ctx,'classes',token)) return;
      const aSamples = biological.filter((row) => row[state.field] === state.first);
      const bSamples = biological.filter((row) => row[state.field] === state.second);
      if (!aSamples.length || !bSamples.length || state.first === state.second) {
        summary.textContent = 'Choose two distinct biological groups with saved raw measurements.';
        results.replaceChildren();
        if (chart.data && window.Plotly) Plotly.purge(chart);
        chart.replaceChildren();
        return;
      }
      const node = branch();
      const mean = (ids,samples) => {
        let total = 0;
        for (const sample of samples) {
          const measurements = rawBySample.get(sampleKey(sample));
          for (const id of ids) {
            const value = number(measurements[id]);
            if (Number.isFinite(value) && value > 0) total += value;
          }
        }
        return total / samples.length;
      };
      const baseA = mean(node.ids,aSamples), baseB = mean(node.ids,bSamples);
      summary.textContent = state.first+' ('+aSamples.length+' samples) vs '+state.second+' ('+bSamples.length+
        ' samples) · '+node.ids.length+' distinct CANOPUS-predicted measured features in this branch.';
      trail.replaceChildren();
      const crumbs = [{name:root.name,depth:0}];
      let ancestor = root;
      for (const [index,name] of state.path.entries()) {
        ancestor = ancestor.children.get(name);
        crumbs.push({name:ancestor.name,depth:index+1});
      }
      for (const crumb of crumbs) {
        if (crumb.depth) trail.append(el('span',{text:'›','aria-hidden':'true'}));
        const button = el('button',{type:'button',class:'text-button',text:crumb.name,
          'aria-current':crumb.depth === state.path.length ? 'page' : 'false'});
        button.addEventListener('click',() => { state.path.length = crumb.depth; update(); });
        trail.append(button);
      }
      const children = [...node.children.values()].map((child) => {
        const a = mean(child.ids,aSamples), b = mean(child.ids,bSamples);
        return {child,a,b,shareA:baseA ? 100*a/baseA : 0,shareB:baseB ? 100*b/baseB : 0};
      });
      const precision = (value) => Number(value.toPrecision(5));
      const level = ['','NPC pathway','NPC superclass','NPC class','Feature ID'][node.level+1] || 'Feature ID';
      const view = table([level,'Features',state.first+' mean raw '+ctx.run.intensity_measure,
        state.second+' mean raw '+ctx.run.intensity_measure,state.second+' / '+state.first,
        state.first+' share (%)',state.second+' share (%)'],
      children.map(({child,a,b,shareA,shareB}) =>
        [child.name,child.ids.length,precision(a),precision(b),a > 0 ? precision(b/a) : 'Not defined',
          precision(shareA),precision(shareB)]),(name) => {
        if (node.level === 3) openFeature(name);
        else { state.path.push(name); update(); }
      },'npc-'+node.level+'-'+ctx.run.hash.slice(0,12)+'.tsv');
      results.replaceChildren(note('Choose a row to drill down. Click a feature ID for its saved evidence and sample intensities.'),view);
      if (chart.data && window.Plotly) Plotly.purge(chart);
      const top = children.sort((left,right) => right.shareA+right.shareB-left.shareA-right.shareB).slice(0,15).reverse();
      const labels = top.map(({child}) => escapePlot(child.name));
      const traces = [[state.first,'shareA'],[state.second,'shareB']].map(([name,key]) =>
        ({type:'bar',orientation:'h',name:escapePlot(name),y:labels,x:top.map((row) => row[key]),
          marker:{color:groupColor(name)},customdata:top.map((row) => row.child.ids.length),
          hovertemplate:'%{y}<br>'+escapePlot(name)+': %{x:.2f}% of parent<br>%{customdata} features<extra></extra>'}));
      chart.style.height = Math.max(390,top.length*32+120)+'px';
      chart.dataset.sampleGroups = 'true';
      plotting(ctx,chart,traces,{barmode:'group',xaxis:{title:'Share of parent mean raw intensity (%)',rangemode:'tozero'},
        yaxis:{automargin:true},margin:{l:Math.min(260,Math.max(110,...labels.map((name) => name.length*7))),r:20,t:25,b:65}});
      if (typeof chart.removeAllListeners === 'function') chart.removeAllListeners('plotly_click');
      chart.on('plotly_click',(event) => {
        const name = top[event.points?.[0]?.pointNumber]?.child.name;
        if (!name) return;
        if (node.level === 3) openFeature(name);
        else { state.path.push(name); update(); }
      });
    }
    fieldSelect.addEventListener('change',() => { state.field = fieldSelect.value; groupOptions(); update(); });
    aSelect.addEventListener('change',() => { state.first = aSelect.value; update(); });
    bSelect.addEventListener('change',() => { state.second = bSelect.value; update(); });
    groupOptions();
    update();
  }

  async function download(ctx,path,button) {
    const original = button.textContent;
    button.disabled = true;
    button.textContent = 'Verifying…';
    try {
      const bytes = await verifiedBytes(ctx,path);
      if (!valid(ctx)) return;
      const type = path.endsWith('.tsv') ? 'text/tab-separated-values' : path.endsWith('.pdf') ? 'application/pdf' : 'image/png';
      const url = URL.createObjectURL(new Blob([bytes],{type}));
      const link = el('a',{href:url,download:path.split('/').pop()});
      document.body.append(link);
      link.click();
      link.remove();
      setTimeout(() => URL.revokeObjectURL(url),60000);
      status('Verified SHA-256 and downloaded '+path+'.');
    } catch (error) {
      if (valid(ctx)) status('Download failed: '+error.message,true);
    } finally {
      if (valid(ctx)) { button.disabled = false; button.textContent = original; }
    }
  }
  function downloadButtons(ctx,paths) {
    if (!paths.length) return note('No corresponding saved files in this run.');
    return el('div',{class:'download-list'},...paths.map((path) => {
      const button = el('button',{type:'button',class:'button secondary small',text:'Download verified '+path.split('.').pop().toUpperCase()});
      button.addEventListener('click',() => download(ctx,path,button));
      return el('div',{class:'download-item'},
        el('div',{},el('span',{class:'mono',text:path}),el('code',{class:'checksum',text:'SHA-256 ' + ctx.run.outputs[path]})),
        button);
    }));
  }
  function renderDownloads(ctx,target) {
    const files = Object.keys(ctx.run.outputs).filter((path) => outputPath.test(path)).sort();
    const tsv = files.filter((path) => path.startsWith('tables/'));
    const plots = files.filter((path) => path.startsWith('plots/'));
    const jobs = catalogue.map((run) =>
      [run.hash,run.created_at,run.dataset_id,run.description ?? 'Not in this catalogue',
        run.recipe_sha256 ?? 'Not in this catalogue',
        run.changes_vs_previous_dataset_run ?? 'Not in this catalogue',run.design_group]);
    target.replaceChildren(
      card('Completed jobs · select a hash to view its saved results',
        table(['Run hash','Completed (UTC)','Dataset ID','Description','Recipe SHA-256','Changes vs previous dataset run','Design group'],
          jobs,(hash) => chooseRun(hash,true), 'completed-jobs-filtered.tsv')),
      card('Selected run',table(['Property','Saved value'],[
        ['Dataset ID',ctx.run.dataset_id],['Run hash',ctx.run.hash],['Description',ctx.run.description ?? 'Not in this catalogue'],
        ['Recipe SHA-256',ctx.run.recipe_sha256 ?? 'Not in this catalogue'],
        ['Changes vs previous dataset run',ctx.run.changes_vs_previous_dataset_run ?? 'Not in this catalogue'],
        ['Created at',ctx.run.created_at],['Schema version',ctx.run.schema_version],
        ['Analysis method','V2 saved statistics'],['Design group',ctx.run.design_group],
        ['Intensity measure',ctx.run.intensity_measure],['Transformation',ctx.run.transformation]],null,'run-properties.tsv')),
      card('Saved TSV tables · SHA-256 verified before download',downloadButtons(ctx,tsv)),
      card('Saved PNG/PDF · original unfiltered plots',downloadButtons(ctx,plots)),
      note('All filtering, dimensions, searches and box plots on this page are descriptive views of published values. Inferential test results are never recomputed.'));
  }
  // Advanced is independent of the selected subset. Only the catalogue's canonical full-batch
  // sealed run supplies its source metadata, through the same SHA-256 verifier as saved views.
  async function renderAdvanced(ctx,target,token) {
    const source = builderSource;
    if (!source) {
      message(target,'A complete original-batch sealed run is not published for this batch. This older catalogue cannot build a new plan; select or publish a canonical full-batch run first.');
      return;
    }
    const run = catalogue.find((item) => item.hash === source.run_hash);
    if (!run || !hex.test(source.run_hash) || source.dataset_id !== run.dataset_id ||
        !Object.hasOwn(run.outputs,'tables/sample_metadata.tsv')) {
      message(target,'The catalogue does not contain the full-batch source metadata and its checksum. A new plan cannot be made safely.',true);
      return;
    }
    if (!ctx.advancedMetadata) {
      const sourceContext = {run,controller:ctx.controller};
      ctx.advancedMetadata = verifiedBytes(sourceContext,'tables/sample_metadata.tsv')
        .then((bytes) => parseTSV(new TextDecoder('utf-8',{fatal:true}).decode(bytes)));
    }
    const full = await ctx.advancedMetadata;
    if (!viewValid(ctx,'advanced',token)) return;
    const roles = source.roles;
    if (!roles?.column || !roles.levels?.sample || !roles.levels?.qc || !roles.levels?.blank ||
        !full.length || !Object.hasOwn(full[0],roles.column) ||
        !/^[A-Za-z0-9_-]+\/[A-Za-z0-9_-]+$/.test(source.dataset_id)) {
      message(target,'The published full-batch source has invalid role or dataset details. No plan can be made safely.',true);
      return;
    }
    const counts = {sample:0,qc:0,blank:0,other:0};
    for (const row of full) {
      const role = Object.entries(roles.levels).find(([,level]) => row[roles.column] === level)?.[0];
      counts[role && Object.hasOwn(counts,role) ? role : 'other']++;
    }
    const samples = full.filter((row) => row[roles.column] === roles.levels.sample);
    const excluded = new Set(['filename','sample_id','sample','analysis_group']);
    const columns = Object.keys(full[0]).filter((column) =>
      column !== roles.column && !excluded.has(column.toLowerCase()) &&
      samples.every((row) => row[column] && row[column].trim() && row[column] !== 'NA') &&
      new Set(samples.map((row) => row[column])).size >= 2 &&
      new Set(samples.map((row) => row[column])).size <= 30);
    if (!columns.length) {
      message(target,'The full-batch sample metadata has no complete categorical biological grouping column.',true);
      return;
    }
    const defaults = source.defaults || {};
    const prep = defaults.preprocessing || {};
    const analyses = defaults.analyses || {};
    if (!ctx.advanced) {
      ctx.advanced = {
        column:columns.includes(source.group_column) ? source.group_column : columns[0],
        groups:new Map(), slug:'merged-groups', description:'Group comparison by '+(columns.includes(source.group_column) ? source.group_column : columns[0]),
        seed:defaults.seed ?? 275,
        preprocessing:{
          blank_filter:{enabled:prep.blank_filter?.enabled ?? true,
            minimum_sample_to_blank_ratio:prep.blank_filter?.minimum_sample_to_blank_ratio ?? 5,
            minimum_blank_prevalence:prep.blank_filter?.minimum_blank_prevalence ?? 0.5},
          qc_rsd_filter:{enabled:prep.qc_rsd_filter?.enabled ?? true,
            threshold_percent:prep.qc_rsd_filter?.threshold_percent ?? 30},
          missing_values:{method:prep.missing_values?.method ?? 'half_minimum'},
          normalization:{method:prep.normalization?.method ?? 'none'},
          transformation:{method:prep.transformation?.method ?? 'log2'},
          scaling:{method:prep.scaling?.method ?? 'pareto'}},
        analyses:{
          pca:{enabled:analyses.pca?.enabled ?? true,components:analyses.pca?.components ?? 5},
          pcoa:{enabled:analyses.pcoa?.enabled ?? true,components:analyses.pcoa?.components ?? 5,
            stage:analyses.pcoa?.stage ?? 'normalized',distance:analyses.pcoa?.distance ?? 'bray'},
          omnibus:{enabled:analyses.omnibus?.enabled ?? true},
          differential:{enabled:analyses.differential?.enabled ?? true},
          plsda:{enabled:analyses.plsda?.enabled ?? false,components:analyses.plsda?.components ?? 3,
            folds:analyses.plsda?.folds ?? 3,repeats:analyses.plsda?.repeats ?? 5,
            permutations:analyses.plsda?.permutations ?? 99}}
      };
    }
    const state = ctx.advanced;
    if (!columns.includes(state.column)) state.column = columns[0];
    const sourceGroups = () => [...new Set(samples.map((row) => row[state.column]))].sort((a,b) => a.localeCompare(b));
    function initializeGroups() {
      state.groups = new Map(sourceGroups().map((name) => [name,{include:true,target:name}]));
    }
    if (!state.groups.size) initializeGroups();
    const choice = (title,value,values,change) => {
      const node = el('select',{},...values.map(([key,label]) => el('option',{value:key,text:label})));
      node.value = value;
      node.addEventListener('change',() => { change(node.value); update(); });
      return el('label',{class:'advanced-field'},el('span',{text:title}),node);
    };
    const numeric = (title,value,min,max,step,change) => {
      const input = el('input',{type:'number',required:'',min:String(min),step:String(step),value:String(value)});
      if (max !== null) input.max = String(max);
      input.addEventListener('input',() => { change(input.value); update(); });
      return el('label',{class:'advanced-field'},el('span',{text:title}),input);
    };
    const toggle = (title,value,change) => {
      const input = el('input',{type:'checkbox'});
      input.checked = value === true;
      input.addEventListener('change',() => { change(input.checked); update(); });
      return el('label',{class:'advanced-toggle'},input,document.createTextNode(title));
    };
    const field = choice('Source group column',state.column,columns.map((name) => [name,name]),(value) => {
      const previousDefault = 'Group comparison by '+state.column;
      state.column = value;
      if (state.description === previousDefault) {
        state.description = 'Group comparison by '+value;
        descriptionInput.value = state.description;
      }
      initializeGroups();
      drawGroups();
    });
    const groupRows = el('tbody');
    const groupPreview = el('div',{class:'advanced-preview','aria-live':'polite'});
    const issues = el('div',{class:'advanced-issues',role:'status','aria-live':'polite'});
    function drawGroups() {
      groupRows.replaceChildren(...sourceGroups().map((name) => {
        const entry = state.groups.get(name);
        const include = el('input',{type:'checkbox','aria-label':'Include '+name});
        include.checked = entry.include;
        include.addEventListener('change',() => { entry.include = include.checked; update(); });
        const targetInput = el('input',{type:'text',value:entry.target,'aria-label':'Merged target label for '+name,required:'',spellcheck:'false'});
        targetInput.addEventListener('input',() => { entry.target = targetInput.value; update(); });
        return el('tr',{},el('th',{scope:'row',text:name}),
          el('td',{text:String(samples.filter((row) => row[state.column] === name).length)}),
          el('td',{},include),el('td',{},targetInput));
      }));
      update();
    }
    const slugInput = el('input',{type:'text',value:state.slug,required:'',pattern:'[a-z][a-z0-9_-]{2,50}',spellcheck:'false'});
    slugInput.addEventListener('input',() => { state.slug = slugInput.value; update(); });
    const slug = el('label',{class:'advanced-field'},el('span',{text:'Output slug (3–51 characters; starts with a letter)'}),slugInput);
    const descriptionInput = el('input',{type:'text',value:state.description,maxlength:'160',
      'aria-label':'Plan description (optional, at most 160 characters)'});
    descriptionInput.addEventListener('input',() => { state.description = descriptionInput.value; update(); });
    const description = el('label',{class:'advanced-field'},
      el('span',{text:'Description (optional, one line, up to 160 characters)'}),descriptionInput);
    const seed = numeric('Random seed',state.seed,0,2147483647,1,(value) => { state.seed = value; });
    const pp = state.preprocessing, an = state.analyses;
    const filterControls = el('div',{class:'advanced-fields'},
      toggle('Blank filter',pp.blank_filter.enabled,(v) => { pp.blank_filter.enabled = v; }),
      numeric('Minimum sample / blank ratio',pp.blank_filter.minimum_sample_to_blank_ratio,0,null,'any',(v) => { pp.blank_filter.minimum_sample_to_blank_ratio = v; }),
      numeric('Minimum blank prevalence (0–1)',pp.blank_filter.minimum_blank_prevalence,0,1,'any',(v) => { pp.blank_filter.minimum_blank_prevalence = v; }),
      toggle('QC RSD filter',pp.qc_rsd_filter.enabled,(v) => { pp.qc_rsd_filter.enabled = v; }),
      numeric('QC RSD threshold (%)',pp.qc_rsd_filter.threshold_percent,0,null,'any',(v) => { pp.qc_rsd_filter.threshold_percent = v; }),
      ...[
        ['missing_values','Missing values',['none','half_minimum']],
        ['normalization','Normalization',['none','total_sum','median','pqn']],
        ['transformation','Transformation',['none','log2','log10','log1p']],
        ['scaling','Scaling',['none','pareto','autoscale']]
      ].map(([key,title,options]) => choice(title,pp[key].method,options.map((name) => [name,name.replaceAll('_',' ')]),(v) => { pp[key].method = v; })));
    const pcoaStage = choice('PCoA matrix stage',an.pcoa.stage,
      ['imputed','normalized','transformed','scaled'].map((v) => [v,v]),(v) => { an.pcoa.stage = v; });
    const pcoaDistance = choice('PCoA distance',an.pcoa.distance,
      ['bray','euclidean','maximum','manhattan','canberra','binary','minkowski'].map((v) => [v,v]),(v) => { an.pcoa.distance = v; });
    const analysisControls = el('div',{class:'advanced-fields'},
      toggle('PCA',an.pca.enabled,(v) => { an.pca.enabled = v; }),
      numeric('PCA components',an.pca.components,1,null,1,(v) => { an.pca.components = v; }),
      toggle('PCoA',an.pcoa.enabled,(v) => { an.pcoa.enabled = v; }),
      numeric('PCoA components',an.pcoa.components,1,null,1,(v) => { an.pcoa.components = v; }),
      pcoaStage,pcoaDistance,
      toggle('Omnibus test',an.omnibus.enabled,(v) => { an.omnibus.enabled = v; }),
      toggle('Differential contrasts (all pairs)',an.differential.enabled,(v) => { an.differential.enabled = v; }),
      toggle('Validated PLS-DA',an.plsda.enabled,(v) => { an.plsda.enabled = v; }),
      numeric('PLS-DA components',an.plsda.components,1,null,1,(v) => { an.plsda.components = v; }),
      numeric('PLS-DA folds',an.plsda.folds,2,null,1,(v) => { an.plsda.folds = v; }),
      numeric('PLS-DA repeats',an.plsda.repeats,2,null,1,(v) => { an.plsda.repeats = v; }),
      numeric('PLS-DA permutations',an.plsda.permutations,19,null,1,(v) => { an.plsda.permutations = v; }));
    const downloadPlan = el('button',{type:'button',class:'button',text:'Download JSON plan'});
    function update() {
      const errors = [];
      const merged = new Map();
      const selected = [];
      for (const name of sourceGroups()) {
        const entry = state.groups.get(name);
        if (!entry.include) continue;
        const target = entry.target.trim();
        if (!/^[A-Za-z][A-Za-z0-9_-]*$/.test(target) || target === 'NA')
          errors.push('Target for '+name+' must start with a letter, use only letters, digits, _ or -, and not be NA.');
        selected.push({source:name,target});
        merged.set(target,(merged.get(target) || 0) + samples.filter((row) => row[state.column] === name).length);
      }
      if (merged.size < 2) errors.push('Select at least two distinct target biological groups.');
      if ([...merged.values()].some((count) => count < 3)) errors.push('Each target biological group needs at least three samples.');
      const smallest = Math.min(...merged.values());
      if (an.plsda.enabled && ([...merged.values()].some((count) => count < 4) || Number(an.plsda.folds) > smallest))
        errors.push('PLS-DA needs at least four samples per group and folds no greater than the smallest group.');
      if (pp.blank_filter.enabled && !counts.blank) errors.push('Blank filtering requires blank samples.');
      if (pp.qc_rsd_filter.enabled && counts.qc < 3) errors.push('QC RSD filtering requires at least three QC samples.');
      if (counts.other) errors.push('The source has '+counts.other+' unrecognized sample roles; review the sealed run before materializing.');
      if (an.pcoa.enabled && an.pcoa.distance === 'bray' && !['imputed','normalized'].includes(an.pcoa.stage))
        errors.push('Bray-Curtis PCoA requires the imputed or normalized nonnegative stage. Select another stage or distance.');
      if (!/^[a-z][a-z0-9_-]{2,50}$/.test(state.slug)) errors.push('Slug must start with a lowercase letter and contain 3–51 lowercase letters, digits, _ or -.');
      if (state.description && (state.description.trim() !== state.description ||
          [...state.description].length > 160 || /[\x00-\x1f\x7f-\x9f]/u.test(state.description)))
        errors.push('Description must be one trimmed line of at most 160 characters.');
      function integer(value,title,min,max = Number.MAX_SAFE_INTEGER) {
        const n = Number(value);
        if (value === '' || !Number.isSafeInteger(n) || n < min || n > max) errors.push(title+' must be an integer from '+min+' to '+max+'.');
        return n;
      }
      function decimal(value,title,max = Infinity) {
        const n = Number(value);
        if (value === '' || !Number.isFinite(n) || n < 0 || n > max) errors.push(title+' must be a number between 0 and '+(max === Infinity ? 'infinity' : max)+'.');
        return n;
      }
      const seedValue = integer(state.seed,'Seed',0,2147483647);
      const blankRatio = decimal(pp.blank_filter.minimum_sample_to_blank_ratio,'Sample / blank ratio');
      const blankPrevalence = decimal(pp.blank_filter.minimum_blank_prevalence,'Blank prevalence',1);
      const qcThreshold = decimal(pp.qc_rsd_filter.threshold_percent,'QC RSD threshold');
      const pcaComponents = integer(an.pca.components,'PCA components',1);
      const pcoaComponents = integer(an.pcoa.components,'PCoA components',1);
      const plsdaComponents = integer(an.plsda.components,'PLS-DA components',1);
      const folds = integer(an.plsda.folds,'PLS-DA folds',2);
      const repeats = integer(an.plsda.repeats,'PLS-DA repeats',2);
      const permutations = integer(an.plsda.permutations,'PLS-DA permutations',19);
      for (const [key,options] of Object.entries({
        missing_values:['none','half_minimum'],normalization:['none','total_sum','median','pqn'],
        transformation:['none','log2','log10','log1p'],scaling:['none','pareto','autoscale']})) {
        if (!options.includes(pp[key].method)) errors.push('Unsupported '+key+' method in source defaults.');
      }
      if (!['imputed','normalized','transformed','scaled'].includes(an.pcoa.stage) ||
          !['bray','euclidean','maximum','manhattan','canberra','binary','minkowski'].includes(an.pcoa.distance))
        errors.push('Unsupported PCoA stage or distance in source defaults.');
      groupPreview.replaceChildren(
        el('p',{text:selected.reduce((sum,item) => sum + samples.filter((row) => row[state.column] === item.source).length,0)+
          ' / '+counts.sample+' biological samples kept · '+counts.qc+' QC and '+counts.blank+' blanks carried through.'}),
        el('ul',{},...[...merged.entries()].sort(([a],[b]) => a.localeCompare(b)).map(([name,count]) =>
          el('li',{text:name+': '+count+' biological samples'}))));
      issues.replaceChildren(...errors.map((error) => note(error,true)));
      downloadPlan.disabled = errors.length > 0;
      return errors.length ? null : {
        schema_version:1,source_run_hash:source.run_hash,slug:state.slug,group_column:state.column,
        ...(state.description ? {description:state.description} : {}),
        groups:selected,seed:seedValue,
        preprocessing:{
          blank_filter:{enabled:pp.blank_filter.enabled,minimum_sample_to_blank_ratio:blankRatio,minimum_blank_prevalence:blankPrevalence},
          qc_rsd_filter:{enabled:pp.qc_rsd_filter.enabled,threshold_percent:qcThreshold},
          missing_values:{method:pp.missing_values.method},normalization:{method:pp.normalization.method},
          transformation:{method:pp.transformation.method},scaling:{method:pp.scaling.method}},
        analyses:{
          pca:{enabled:an.pca.enabled,components:pcaComponents},
          pcoa:{enabled:an.pcoa.enabled,components:pcoaComponents,stage:an.pcoa.stage,distance:an.pcoa.distance},
          omnibus:{enabled:an.omnibus.enabled},differential:{enabled:an.differential.enabled},
          plsda:{enabled:an.plsda.enabled,components:plsdaComponents,folds,repeats,permutations}},
        contrasts:'all'};
    }
    downloadPlan.addEventListener('click',() => {
      const plan = update();
      if (!plan) return;
      const blob = new Blob([JSON.stringify(plan,null,2)+'\n'],{type:'application/json'});
      const url = URL.createObjectURL(blob);
      const link = el('a',{href:url,download:plan.slug+'-plan.json'});
      document.body.append(link);
      link.click();
      link.remove();
      setTimeout(() => URL.revokeObjectURL(url),60000);
    });
    const command = el('pre',{class:'advanced-commands'});
    function instructions() {
      if (!/^[a-z][a-z0-9_-]{2,50}$/.test(state.slug)) {
        command.textContent = 'Enter a valid slug to display local commands.';
        return;
      }
      const batch = '../../didier-reinhardt-group/docs/'+source.dataset_id;
      const base = batch+'/recipes/'+state.slug;
      command.textContent = '# From biostat_toolbox/v2; replace the --plan path with your downloaded file\n'+
        'rig run -r 4.6.1 -f ../scripts/materialize_stats_recipe.R -- --batch '+batch+
        ' --plan /path/to/'+state.slug+'-plan.json\n\n'+
        '# After materialization (it prints these commands too):\n'+
        'rig run -r 4.6.1 -f run.R -- validate --dataset '+base+'/dataset.yaml --recipe '+base+'/recipe.yaml\n'+
        'rig run -r 4.6.1 -f run.R -- plan --dataset '+base+'/dataset.yaml --recipe '+base+'/recipe.yaml\n'+
        'rig run -r 4.6.1 -f run.R -- run --dataset '+base+'/dataset.yaml --recipe '+base+'/recipe.yaml';
    }
    slugInput.addEventListener('input',instructions);
    target.replaceChildren(
      card('Complete original batch (not the selected run)',
        note(source.dataset_id+' · source run '+source.run_hash+' · SHA-256 verified tables/sample_metadata.tsv'),
        metrics([['All rows',full.length],['Biological samples',counts.sample],['QC',counts.qc],['Blanks',counts.blank],['Other roles',counts.other]]),
        note('The current run picker affects saved views only. Advanced always starts with the canonical original metadata; candidate annotations are not sample metadata.')),
      card('1 · Choose and merge biological groups',field,
        el('div',{class:'table-wrap'},el('table',{class:'advanced-group-table'},
          el('thead',{},el('tr',{},...['Original group','Samples','Include','Merged target label'].map((heading) => el('th',{scope:'col',text:heading})))),
          groupRows)),
        el('h4',{text:'Preview: resulting groups'}),groupPreview,issues),
      card('2 · Name the plan',el('div',{class:'advanced-fields'},slug,description,seed)),
      card('3 · V2 preprocessing',note('The original QC and blank rows remain available for enabled filters.'),filterControls),
      card('4 · V2 analyses',analysisControls),
      card('5 · Download and run locally',downloadPlan,
        note('This is a configuration plan only. Save the JSON file locally, then run the trusted repository materializer from v2/. It verifies original input SHA-256, validates the generated dataset and never overwrites an existing output directory. No browser-side analysis or server write occurs.'),
        command));
    drawGroups();
    instructions();
  }
  if (detachedMode) {
    document.body.classList.add('detached-view');
    $('detached-toolbar').hidden = false;
    $('detached-attach').disabled = !window.opener;
    for (const tab of tabs) $(tab).hidden = tab !== 'features';
  }
  boot();
})();
