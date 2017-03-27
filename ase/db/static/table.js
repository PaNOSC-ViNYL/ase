
//
// On load
//

function BodyOnLoad(id, query)
{
    //console.log(id);

    if(sessionStorage.getItem("cid") !== id.toString())
    {
        // create and store the session id
        sessionStorage.cid = id;

        // reset collapse if a new session id is given
        sessionStorage.removeItem('collapseExtraSearch');
      
        // 
        // setup suggestions
        SetupSuggestions();
      
        //
        // setup controls
        CopyCtrl();
    }

    ns.Init(query);
    
    document.getElementById("formula-result").focus();
}

//
// Requests
//

function SetupSuggestions()
{
    var xhr = new XMLHttpRequest(),
    method = "GET",
    url = "/formulas";

    xhr.open(method, url, false);
    xhr.onload = function() 
    {
        if(xhr.readyState === XMLHttpRequest.DONE && xhr.status === 200) 
        {
            //console.log(xhr.responseText);
            sessionStorage.setItem("formula", xhr.responseText);
        }
    };

    xhr.send();
}

function CopyCtrl()
{
    var xhr = new XMLHttpRequest(),
    method = "GET",
    url = "/special_keys";

    xhr.open(method, url, false);
    xhr.onload = function() 
    {
        if(xhr.readyState === XMLHttpRequest.DONE && xhr.status === 200) 
        {
            var ctrl = JSON.parse(xhr.responseText);
            //console.log(ctrl);

            var cList = [];
            var cType = [];

            for(var iter in ctrl)
            {
                cList.push(ctrl[iter][0]);
                cType.push(ctrl[iter][1]);
            }

            sessionStorage.setItem("ctrlKeys", JSON.stringify(cList));
            sessionStorage.setItem("ctrlType", JSON.stringify(cType));
        }
    };

    xhr.send();   
}

//
// remove multiple instances from array
//

function uniq(a)
{
    return a.sort().filter(function(item, pos, ary) 
    {
        return !pos || item != ary[pos - 1];
    })
}

//
// input controls calls this function on update
//

function ControlFunction(type, element, value)
{
    if(value === 'Reset')
        ns.SetField(element, '');
    else
        ns.SetField(element, value);
    
    document.getElementById("formula-result").focus(); 
}
        
//
// document ready
//

$(document).ready(function()
{
    //
    // Setup auto complete
    nsAuto.Init();

    //
    // simple search bar (ControlFunction)
    document.getElementById("formula-result").onchange = function()
    {
        ns.SetField('formula', this.value);
    };

    //
    // autocomplete
    //var autoSuggestions = JSON.parse(sessionStorage.getItem("formula"));
    //console.log(autoSuggestions);
    $( "#formula-result" ).autocomplete(
    {
        source: nsAuto.Get()
    });

    //
    // tooltips
    $('[data-toggle="tooltip"]').tooltip();

    //
    // collapse well storage
    $(".collapse").on('shown.bs.collapse', function()
    {
        if (this.id) 
        {
            sessionStorage[this.id] = 'show';
        }
    }).on('hidden.bs.collapse', function()
    {
        if (this.id) 
        {
            sessionStorage.removeItem(this.id);
        }    
    }).each(function()
    {
        if(this.id && sessionStorage[this.id] === 'show')
        {
            $(this).collapse('show');
        }
    });
})

//
// namespace for autosuggestions
//

var nsAuto = (function()
{
    var m_autoSuggestions;

    function Init()
    {
        this.m_autoSuggestions = JSON.parse(sessionStorage.getItem("formula"));
    }

    function Get()
    {
        return this.m_autoSuggestions;
    }

    //
    // public function interface
    return {
        Init : Init,
        Get : Get
    };
})()
//
// class for input control
//

class inputCtrl
{
    constructor(key, type, index)
    {
        this.m_key = key;
        this.m_type = type; // ['text', 'COMBO', 'INTERVAL', 'CHECK']
        this.m_access = -1;

        if(this.m_type === "COMBO" || this.m_type === "text")
            this.m_value = '';
        else if(this.m_type === "CHECK")
            this.m_value = 0;
        else if(this.m_type === "INTERVAL")
            this.m_value = [null, null];
    }

    GetQueryString()
    {
        if(this.m_type === "text")
        {
            var formulas = nsAuto.Get();

            if(formulas.indexOf(this.m_value) === -1)
            {
                return this.m_value;
            }
            else
            {
                return this.m_key + '=' + this.m_value;
            }
        }
        if(this.m_type === "COMBO")
        {
            return this.m_key + '=' + this.m_value;
        }
        else if(this.m_type === "CHECK")
        {
            if(this.m_value === 1)
                return this.m_key;
            else
                return '';
        }
        else if(this.m_type === "INTERVAL")
        {
            var str = '';
            if(this.m_value[0] !== null)
                str += this.m_value[0] + '<=';

            str += this.m_key;

            if(this.m_value[1] !== null)
                str += '<=' + this.m_value[1];

            return str;
        }
        
        return '';
    }

    IsValid()
    {
        if(this.m_type === "COMBO" || this.m_type === "text")
        {
            if(this.m_value == '')
                return false;
            else
                return true;
        }
        else if(this.m_type === "CHECK")
        {
            return this.m_value;
        }
        else if(this.m_type === "INTERVAL")
        {
            if(this.m_value[0] == null && this.m_value[1] == null)
                return false;
            return true;
        }
    }

    ForceAssign(tokenID, token)
    {
        // the ForceAssign function should not be abused

        if(this.m_type === "text")
        {
            this.m_access = tokenID;
            this.m_value = token;
            return true;
        }

        return false;
    }

    TryAssign(tokenID, token)
    {        
        if(this.m_type === "COMBO" || this.m_type === "text")
        {
            // we want a perfect match here
            var checkstr = token.substr(0, this.m_key.length + 1);
            				
			if(checkstr === this.m_key + '=')
            {
                this.m_access = tokenID;
                this.m_value = token.substr(this.m_key.length + 1, token.length);

                return true;
            }        
        }
        else if(this.m_type === "CHECK")
        {
            if(this.m_key === token)
            {
                this.m_access = tokenID;
                this.m_value = 1;
                return true;
            }
        }
        else if(this.m_type === "INTERVAL")
        {
            // the key is inside the token
            var keyPos = ns.getIndicesOf(this.m_key, token);
            var indices = ns.getIndicesOf('<=', token);

            if(indices.length === 1)
            {
                if(indices[0] < keyPos[0])
                {
                    // check layout
                    var chlayout = ns.getIndicesOf('<=' + this.m_key, token);
                
                    if(chlayout.length === 1)
                    {
                        var v1 = token.substr(0, indices[0]);

                        this.m_access = tokenID;
                        this.m_value = [v1, null];

                        return true;
                    }
                }
                else
                {
                    // check layout
                    var chlayout = ns.getIndicesOf(this.m_key + '<=', token);
                
                    if(chlayout.length === 1)
                    {
                        var pos = indices[0]+2;
                        var v2 = token.substr(pos, token.length-pos);

                        this.m_access = tokenID;
                        this.m_value = [null, v2];

                        return true;
                    }
                }
            
                return false;
            }
            else if(indices.length === 2)
            {
                // check layout
                var chlayout = ns.getIndicesOf('<=' + this.m_key + '<=', token);
                
                if(chlayout.length === 1)
                {
                    // extract values
                    var v1 = token.substr(0, indices[0]);
                    var pos = indices[1]+2;
                    var v2 = token.substr(pos, token.length-pos);

                    this.m_access = tokenID;
                    this.m_value = [v1, v2];

                    return true;
                }

                return false;
            }
        }

        return false;
    }

    SetOutput()
    {
        var element = this.m_key + '-result';

        if(this.m_type === "text")
			document.getElementById(element).value = this.m_value;
        else if(this.m_type === "CHECK")
			document.getElementById(element).checked = this.m_value;
        else if(this.m_type === "INTERVAL")
        {
            if(this.m_value[0] !== null)
			    document.getElementById(this.m_key + '-l').value = this.m_value[0];
            if(this.m_value[1] !== null)
			    document.getElementById(this.m_key + '-r').value = this.m_value[1];
        }
		else
			document.getElementById(element).innerHTML = this.m_value;
    }
};

//
// ns holds the functionality for the generation of the advanced search string
//

var ns = (function()
{
    //
	// local variables in namespace
    //

    // search string
	var m_query;

    // tokens
	var m_qlist;
	var m_recognized;

    // controls
    var m_control = [];

	//
	// public functions
    function Init(query)
	{
		m_query = query;

		// set query in forms
		document.getElementById('searchstr').value = query;

        //
        // setup supported controls
        SetupControls();

        //
        // retrive tokens from search string
        ExtractTokens();

		//
        // init controls as specified in the search string
		SetLoaded();

        //
        // write query
        UpdateSearchQuery();
	}

    function GetControlID(key)
    {
        for(i=0; i<m_control.length; ++i)
        {
            if(m_control[i].m_key == key)
                return i;
        }
        return -1;
    }
    
    function ClearVariable(index)
    {
        m_qlist.splice(m_control[index].m_access, 1);
		m_recognized.splice(m_control[index].m_access, 1);

		for(i=0;i<m_control.length;++i)
		{
			if(m_control[i].m_access > m_control[index].m_access)
				m_control[i].m_access -= 1;
		}

		m_control[index].m_access = -1;
    }

    function CreateVariable(index, key, value)
    {
        m_qlist.push(m_control[index].GetQueryString());
        m_recognized.push(index);

		m_control[index].m_value = value;
		m_control[index].m_access = m_qlist.length-1;
    }

	function SetField(key, value)
	{
        var index = GetControlID(key);

        if(index === -1)
            return;

        // update the value of th field
        m_control[index].m_value = value;

		if(m_control[index].m_access !== -1)
		{
			// if the value is empty we remove the variable
            if(m_control[index].IsValid() == false)
                ClearVariable(index);    
		}
		else
		{
			// field has not been set, create it
            if(m_control[index].IsValid() == true)
                CreateVariable(index, key, value);
		}

		// update GUI if it is not a textfield
        m_control[index].SetOutput();

		UpdateSearchQuery();
	}

	//
	// private functions (stays invisible)

    function SetupControls()
    {
        var key = JSON.parse(sessionStorage.getItem("ctrlKeys"));
        var type = JSON.parse(sessionStorage.getItem("ctrlType"));

        //console.log(key);
        //console.log(type);

        // allways add the simple searchfield
        m_control.push(new inputCtrl('formula', 'text', key.length));

        for(i=0; i<key.length; i++)
        {
            m_control.push(new inputCtrl(key[i], type[i], i));
        }

        //console.log(m_control);
    }

    function ExtractTokens()
    {
        m_qlist = [];
		m_recognized = [];

		var seperationIndices = getIndicesOf(",", m_query);
		seperationIndices.unshift(-1);
		seperationIndices.push(m_query.length);
			
		for(i=1; i<seperationIndices.length; i++)
		{
			var ss = m_query.substr(seperationIndices[i-1]+1,seperationIndices[i]-seperationIndices[i-1]-1);

            ss = ss.trim();

			if(ss !== '')
			{
				m_qlist.push(ss);
				m_recognized.push(-1);
			}
		}

        //console.log(m_qlist);
    }

    function SetLoaded()
	{
        //
        // For each control find the first appearence in qlist that matches it and try to assign it

		for(j=0; j<m_control.length; j++)
		{
			for(i=0; i<m_qlist.length; i++)
			{
				var indices = getIndicesOf(m_control[j].m_key, m_qlist[i]);                
				
				if(indices.length === 1)
				{
                    if(m_control[j].TryAssign(i, m_qlist[i]) == true)
                    {
						m_recognized[i] = j;
						break;
					}
				}
			}
		}

        // 
        // if the simple search field is not assigned to a token
        // we collect all the unhandled/unrecognized tokens
        // and assign them

        if(m_control[0].m_access === -1)
        {
            var srhQ = "";
            for(var i=0; i<m_qlist.length; i++)
            {
                if(m_recognized[i] === -1)
                    srhQ += m_qlist[i] + ",";
            }

            //
            // remove the unrecognized tokens
            for (var i = m_qlist.length-1; i >= 0; i--)
            {
                if(m_recognized[i] === -1)
                {
                    m_qlist.splice(i, 1);
                    m_recognized.splice(i, 1);
                }
            }

            //
            // add the token if it is non-empty
            if(srhQ !== "")
            {
                srhQ = srhQ.substr(0, srhQ.length-1);
                var index = m_qlist.length;
                m_qlist.push(srhQ);
                m_recognized.push(0);
                m_control[0].ForceAssign(index, m_qlist[index]);                
            }
        }

        //
        // update output controls on the webpage

		for(i=0; i<m_control.length; ++i)
		{
            if(m_control[i].m_access !== -1)
            {
			    m_control[i].SetOutput();
            }
		}
	}

	function UpdateSearchQuery()
	{
		var query = '';

		for(i=0; i<m_qlist.length; i++)
		{
			if(m_recognized[i] === -1)
			{
				query += m_qlist[i];
			}
			else
			{
				if(m_control[m_recognized[i]].m_value === '')
				{
					continue;
				}
				query += m_control[m_recognized[i]].GetQueryString();
			}

			if(i !== m_qlist.length-1)
			{
				query += ',';
			}
		}

		document.getElementById('searchstr').value = query;
	}

    function getIndicesOf(searchStr, str, caseSensitive) 
	{
		var searchStrLen = searchStr.length;
		if (searchStrLen == 0) {
			return [];
		}
		var startIndex = 0, index, indices = [];
		if (!caseSensitive) {
			str = str.toLowerCase();
			searchStr = searchStr.toLowerCase();
		}
		while ((index = str.indexOf(searchStr, startIndex)) > -1) {
			indices.push(index);
			startIndex = index + searchStrLen;
		}
		return indices;
	}

    //
    //
    // stuff visible in namespace object
    return {
        Init : Init,
        getIndicesOf : getIndicesOf,
        SetField : SetField
    };
})()